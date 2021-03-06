#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  ---------------------------------------------------------------------
#
#  _____    _      _              _         _____ _____
# | ____|__| | ___| |_      _____(_)___ ___|  ___| ____|
# |  _| / _` |/ _ \ \ \ /\ / / _ \ / __/ __| |_  |  _|
# | |__| (_| |  __/ |\ V  V /  __/ \__ \__ \  _| | |___
# |_____\__,_|\___|_| \_/\_/ \___|_|___/___/_|   |_____|
#
#
#  Unit of Strength of Materials and Structural Analysis
#  University of Innsbruck,
#  2017 - today
#
#  Matthias Neuner matthias.neuner@uibk.ac.at
#
#  This file is part of EdelweissFE.
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  The full text of the license can be found in the file LICENSE.md at
#  the top level directory of EdelweissFE.
#  ---------------------------------------------------------------------
"""
Created on Sun Jan 15 14:22:48 2017

@author: Matthias Neuner

Output manager for Ensight exports.
If loaded, it automatically exports all elSets as Ensight parts.
For each part, perNode and perElement results can be exported, which are imported from fieldOutputs.

ATTENTION: 
    fieldOutputs for perNode results must be defined on elSets instead of a nodeSet.

Datalines:
"""

documentation = {
    "configuration": "overwrite=y|n|yes|no|True|False|",
    "create": "perNode|perElement",
    "fieldOutput": "name of the result, defined on an elSet (also for perNode results)",
    "name": "(optional), standard = fieldOutputs name",
    "intermediateSaveInterval": 'step option in category "Ensight": save .case file every N increments',
}

from fe.outputmanagers.outputmanagerbase import OutputManagerBase

import os
import datetime
import numpy as np
from collections import defaultdict, OrderedDict
from distutils.util import strtobool
from fe.utils.misc import stringDict
from fe.utils.meshtools import disassembleElsetToEnsightShapes
import fe.config.phenomena
from fe.utils.math import evalModelAccessibleExpression


def writeCFloat(f, ndarray):
    np.asarray(ndarray, dtype=np.float32).tofile(f)


def writeCInt(f, ndarray):
    np.asarray(ndarray, dtype=np.int32).tofile(f)


def writeC80(f, string):
    np.asarray(string, dtype="a80").tofile(f)


variableTypes = {"scalar": 1, "vector": 3, "rotation vector": 3, "tensor": 9}

ensightPerNodeVariableTypes = {1: "scalar per node", 3: "vector per node", 6: "tensor per node", 9: "tensor9 per node"}

ensightPerElementVariableTypes = {
    1: "scalar per element",
    3: "vector per element",
    6: "tensor per element",
    9: "tensor asym per element",
}


class EnsightUnstructuredPart:
    """define an unstructured part, by a list of nodes and a dictionary of elements.
    Each dictionary entry consists of a list of tuples of elementlabel and nodelist:
    {strElementType : [ ( intLabel = None, [nodeList] ) ]}"""

    def __init__(
        self,
        description,
        partNumber,
        nodes,
        elementTree,
    ):

        self.structureType = "coordinates"
        self.elementTree = elementTree

        self.nodes = nodes
        self.nodeLabels = np.asarray([node.label for node in nodes], np.int32)

        self.description = description  # string, describing the part; max. 80 characters
        self.partNumber = partNumber
        self.nodeCoordinateArray = np.asarray([node.coordinates for node in nodes])

    def writeToFile(self, binaryFileHandle, printNodeLabels=True, printElementLabels=True):

        if len(self.nodeCoordinateArray.shape) > 1 and self.nodeCoordinateArray.shape[1] < 3:
            extendTo3D = True
        else:
            extendTo3D = False

        nNodes = self.nodeCoordinateArray.shape[0]
        f = binaryFileHandle

        writeC80(f, "part")
        writeCInt(f, self.partNumber)
        writeC80(f, self.description)
        writeC80(f, "coordinates")
        writeCInt(f, nNodes)

        if printNodeLabels and self.nodeLabels is not None:
            writeCInt(f, self.nodeLabels)

        writeCFloat(f, self.nodeCoordinateArray.T)

        if extendTo3D:
            writeCFloat(f, np.zeros(self.nodeCoordinateArray.shape[0] * (3 - self.nodeCoordinateArray.shape[1])))

        for elemType, elements in self.elementTree.items():
            writeC80(f, elemType)
            writeCInt(f, len(elements))
            if printElementLabels:
                writeCInt(f, np.asarray([element.elNumber for element in elements.keys()], np.int32))

            for nodeIndices in elements.values():
                writeCInt(f, np.asarray(nodeIndices, np.int32)[:] + 1)


# class EnsightPointPart(EnsightUnstructuredPart):
#    """derived UnstructuredPart to represent a single point"""
#    def __init__(self, partNumber, coordinates, description="single_point"):
#        super().__init__(description, partNumber, {"point" : [(1, [0])] }, [ {'label': 1, 'coords':coordinates, } ])


class EnsightTimeSet:
    """defines a set which may be used by EnsightGeometry, EnsightStructuredPart, EnsightUnstructuredPart and is written into the case file"""

    def __init__(
        self, number=1, description="timeStepDesc", fileNameStartNumber=0, fileNameNumberIncrement=1, timeValues=None
    ):
        self.number = number
        self.description = description
        self.fileNameStartNumber = fileNameStartNumber
        self.fileNameNumberIncrement = fileNameNumberIncrement
        self.timeValues = timeValues if timeValues is not None else []


class EnsightGeometryTrend:
    """container class for the time dependent evolution of the geometry,
    establishes the connection between the geometry entities and a EnsightTimeSet"""

    def __init__(self, ensightTimeSet, ensightGeometryList=None):
        self.timeSet = ensightTimeSet
        self.geometryList = ensightGeometryList if ensightGeometryList is not None else []


class EnsightGeometry:
    """container class for one or more EnsightParts at a certain time state, handles also the file writing operation"""

    def __init__(
        self,
        name="geometry",
        descriptionLine1="",
        descriptionLine2="",
        ensightPartList=None,
        nodeIdOption="given",
        elementIdOption="given",
    ):
        self.name = name
        self.descLine1 = descriptionLine1
        self.descLine2 = descriptionLine2
        self.partList = ensightPartList if ensightPartList is not None else []
        self.nodeIdOption = nodeIdOption
        self.elementIdOption = elementIdOption

    def writeToFile(self, fileHandle):
        f = fileHandle
        writeC80(f, self.descLine1)
        writeC80(f, self.descLine2)
        writeC80(f, "node id " + self.nodeIdOption)
        writeC80(f, "element id " + self.elementIdOption)

        if self.nodeIdOption == "given" or self.nodeIdOption == "ignore":
            printNodeLabels = True
        else:  # assign or off
            printNodeLabels = False

        if self.elementIdOption == "given" or self.nodeIdOption == "ignore":
            printElementLabels = True
        else:  # assign or off
            printElementLabels = False

        for part in self.partList:
            part.writeToFile(f, printNodeLabels, printElementLabels)


class EnsightVariableTrend:
    """container class for the time dependent evolution of one variable,
    establishes the connection between EnsightVariable entities and a EnsighTimeSet"""

    def __init__(
        self,
        ensightTimeSet,
        variableName,
        ensightVariableList=None,
        variableType="scalar per node",
        description="variableTrendDescription",
    ):
        self.timeSet = ensightTimeSet
        self.variableName = variableName
        self.variableList = ensightVariableList if ensightVariableList is not None else []
        self.variableType = variableType
        self.description = description


class EnsightPerNodeVariable:
    """container class for data for one certain variable, defined for one or more parts (classification by partID), at a certain time state.
    For each part the structuretype ("coordinate" or "block") has to be defined.
    Each part-variable assignment is defined by a dictionary entry of type: { EnsightPart: np.array(variableValues) }"""

    def __init__(self, name, variableDimension, ensightPartsDict=None):
        self.name = name
        self.description = name
        self.partsDict = ensightPartsDict or {}  # { EnsightPart: np.array(variableValues) }
        self.variableDimension = variableDimension
        self.varType = ensightPerNodeVariableTypes[variableDimension]

    def writeToFile(
        self,
        fileHandle,
    ):
        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, (structureType, values) in self.partsDict.items():
            writeC80(f, "part")
            writeCInt(f, ensightPartID)
            writeC80(f, structureType)
            writeCFloat(f, values.T)
            if values.shape[1] < self.variableDimension:
                writeCFloat(f, np.zeros((values.shape[0], self.variableDimension - values.shape[1])))


class EnsightPerElementVariable:
    """Container class for data for one certain variable, defined for one or more parts (classification by partID), at a certain time state.
    For each part the structuretype ("coordinate" or "block") has to be defined.
    Each part-variable assignment is defined by a dictionary entry of type: { EnsightPart: np.array(variableValues) }"""

    def __init__(
        self,
        name,
        variableDimension,
        ensightPartsDict=None,
    ):
        self.name = name
        self.description = name
        self.partsDict = ensightPartsDict or {}
        self.varType = ensightPerElementVariableTypes[variableDimension]

    def writeToFile(self, fileHandle):

        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, elTypeDict in self.partsDict.items():
            writeC80(f, "part")
            writeCInt(f, ensightPartID)
            for elType, values in elTypeDict.items():
                writeC80(f, elType)
                writeCFloat(f, values.T)


class EnsightChunkWiseCase:
    def __init__(self, caseName, directory="", writeTransientSingleFiles=True):
        self.directory = directory
        self.caseName = caseName
        self.caseFileNamePrefix = os.path.join(directory, caseName)
        self.writeTransientSingleFiles = writeTransientSingleFiles
        self.timeAndFileSets = {}
        self.geometryTrends = {}
        self.variableTrends = {}
        self.fileHandles = {}

        if not os.path.exists(self.caseFileNamePrefix):
            os.mkdir(self.caseFileNamePrefix)

    def setCurrentTime(self, timeAndFileSetNumber, timeValue):
        if not timeAndFileSetNumber in self.timeAndFileSets:
            self.timeAndFileSets[timeAndFileSetNumber] = EnsightTimeSet(timeAndFileSetNumber, "noDesc", 0, 1)
        tfSet = self.timeAndFileSets[timeAndFileSetNumber]
        tfSet.timeValues.append(timeValue)

    def writeGeometryTrendChunk(self, ensightGeometry, timeAndFileSetNumber=1):

        if ensightGeometry.name not in self.fileHandles:
            fileName = os.path.join(
                self.caseFileNamePrefix,
                ensightGeometry.name + ".geo",
            )

            self.fileHandles[ensightGeometry.name] = open(fileName, mode="wb")

        f = self.fileHandles[ensightGeometry.name]

        if not ensightGeometry.name in self.geometryTrends:
            self.geometryTrends[ensightGeometry.name] = timeAndFileSetNumber
            writeC80(f, "C Binary")

        if self.writeTransientSingleFiles:
            writeC80(f, "BEGIN TIME STEP")
            ensightGeometry.writeToFile(f)
            writeC80(f, "END TIME STEP")

    def writeVariableTrendChunk(self, ensightVariable, timeAndFileSetNumber=2):

        if ensightVariable.name not in self.fileHandles:

            fileName = os.path.join(self.caseFileNamePrefix, ensightVariable.name + ".var")

            self.fileHandles[ensightVariable.name] = open(fileName, mode="wb")

        f = self.fileHandles[ensightVariable.name]

        if not ensightVariable.name in self.variableTrends:
            self.variableTrends[ensightVariable.name] = timeAndFileSetNumber, ensightVariable.varType
            writeC80(f, "C Binary")

        if self.writeTransientSingleFiles:
            writeC80(f, "BEGIN TIME STEP")
            ensightVariable.writeToFile(f)
            writeC80(f, "END TIME STEP")

    def finalize(self, replaceTimeValuesByEnumeration=True, closeFileHandes=True):
        caseFName = self.caseFileNamePrefix + ".case"

        if closeFileHandes:
            for f in self.fileHandles.values():
                f.close()

        with open(caseFName, mode="w") as cf:
            cf.write("FORMAT\n")
            cf.write("type: ensight gold\n")

            cf.write("TIME\n")
            for setNum, timeSet in self.timeAndFileSets.items():
                cf.write("time set: " + str(setNum) + " noDesc\n")
                cf.write("number of steps: " + str(len(timeSet.timeValues)) + "\n")
                cf.write("filename start number: " + str(timeSet.fileNameStartNumber) + "\n")
                cf.write("filename increment: " + str(timeSet.fileNameNumberIncrement) + "\n")
                cf.write("time values: ")
                for i, timeVal in enumerate(timeSet.timeValues):
                    if not replaceTimeValuesByEnumeration:
                        cf.write("{:1.8e}".format(timeVal) + "\n")
                    else:
                        cf.write("{:}".format(i) + "\n")

            if self.writeTransientSingleFiles:
                cf.write("FILE\n")
                for timeSet in self.timeAndFileSets.values():
                    cf.write("file set: {:} \n".format(timeSet.number))
                    cf.write("number of steps: {:} \n".format(len(timeSet.timeValues)))

            cf.write("GEOMETRY\n")
            for geometryName, tAndFSetNum in self.geometryTrends.items():
                cf.write("model: {:} \n".format(os.path.join(self.caseFileNamePrefix, geometryName + ".geo")))

            cf.write("VARIABLE\n")
            for variableName, (tAndFSetNum, variableType) in self.variableTrends.items():
                cf.write(
                    "{:}: {:} {:} {:} {:}.var\n".format(
                        variableType,
                        tAndFSetNum,
                        tAndFSetNum,
                        variableName,
                        os.path.join(self.caseFileNamePrefix, variableName),
                    )
                )


def createUnstructuredPartFromElementSet(setName, elementSet, partID):
    """Determines the element and node list for an Ensightpart from an
    element set. The reduced, unique node set is generated, as well as
    the element to node index mapping for the ensight part."""

    nodeCounter = 0
    partNodes = OrderedDict()  # node -> index in nodelist
    elementDict = defaultdict(OrderedDict)
    for element in elementSet:
        elNodeIndices = []
        for node in element.nodes:
            # if the node is already in the dict, get its index,
            # else insert it, and get the current idx = counter. increase the counter
            idx = partNodes.setdefault(node, nodeCounter)
            elNodeIndices.append(idx)
            if idx == nodeCounter:
                # the node was just inserted, so increase the counter of inserted nodes
                nodeCounter += 1
        elementDict[element.ensightType][element] = elNodeIndices

    return EnsightUnstructuredPart(setName, partID, partNodes.keys(), elementDict)


class OutputManager(OutputManagerBase):
    identification = "Ensight Export"

    def __init__(self, name, definitionLines, jobInfo, modelInfo, fieldOutputController, journal, plotter):
        self.name = name

        self.timeAtLastOutput = -1e16
        self.minDTForOutput = -1e16
        self.finishedSteps = 0
        self.intermediateSaveInterval = 10
        self.intermediateSaveIntervalCounter = 0
        self.domainSize = modelInfo["domainSize"]
        self.fieldOutputController = fieldOutputController
        self.journal = journal

        self.transientTAndFSetNumber = 1
        self.staticTAndFSetNumber = 2

        exportName = "{:}_{:}".format(name, datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))

        for defLine in definitionLines:
            definition = stringDict(defLine)
            if "configuration" in definition:
                if strtobool(definition.get("overwrite", "True")):
                    exportName = "{:}".format(name)

        self.elSetToEnsightPartMappings = {}
        self.ensightCase = EnsightChunkWiseCase(exportName)

        self.ensightCase.setCurrentTime(self.staticTAndFSetNumber, 0.0)

        self.transientPerNodeJobs = []
        self.transientPerElementJobs = []

        elementSets = modelInfo["elementSets"]

        elSetParts = []
        partCounter = 1
        for setName, elSet in elementSets.items():
            elSetPart = createUnstructuredPartFromElementSet(setName, elSet, partCounter)
            self.elSetToEnsightPartMappings[setName] = elSetPart
            elSetParts.append(elSetPart)
            partCounter += 1

        geometry = EnsightGeometry("geometry", "Edelweiss_FE", "*export*", ensightPartList=elSetParts)

        geometryTimesetNumber = None
        self.ensightCase.writeGeometryTrendChunk(geometry, geometryTimesetNumber)

        for defLine in definitionLines:
            definition = stringDict(defLine)

            # standard, transient jobs accessing the fieldoutput:
            if "fieldOutput" in definition:
                varType = definition["create"]

                fieldOutput = fieldOutputController.fieldOutputs[definition["fieldOutput"]]
                if fieldOutput.domainType != "elSet":
                    raise Exception("Ensight output can only operate on fieldOutputs defined on elSets!")

                if varType == "perNode":
                    perNodeJob = {}
                    perNodeJob["fieldOutput"] = fieldOutput
                    perNodeJob["name"] = definition.get("name", perNodeJob["fieldOutput"].name).replace(" ", "_")
                    perNodeJob["part"] = self.elSetToEnsightPartMappings[perNodeJob["fieldOutput"].nSetName]
                    field = perNodeJob["fieldOutput"].field
                    perNodeJob["varSize"] = variableTypes[fe.config.phenomena.phenomena[field]]
                    perNodeJob["dimensions"] = fe.config.phenomena.getFieldSize(field, self.domainSize)
                    self.transientPerNodeJobs.append(perNodeJob)

                if varType == "perElement":
                    perElementJob = {}
                    perElementJob["fieldOutput"] = fieldOutput
                    perElementJob["part"] = self.elSetToEnsightPartMappings[perElementJob["fieldOutput"].elSetName]

                    # TODO: don't do this for each output job, but only once!
                    perElementJob["elementsOfShape"] = disassembleElsetToEnsightShapes(
                        perElementJob["fieldOutput"].elSet
                    )
                    perElementJob["name"] = definition.get("name", perElementJob["fieldOutput"].name).replace(" ", "_")
                    self.transientPerElementJobs.append(perElementJob)

            # one time, static jobs, which may output model data; they are executed immediately
            elif "modeldata" in definition:
                varType = definition["create"]
                if varType == "perElement":
                    perElementJob = {}
                    name = definition["name"]
                    elSet = modelInfo["elementSets"][definition["elSet"]]
                    part = self.elSetToEnsightPartMappings[definition["elSet"]]

                    # TODO: don't do this for each output job, but only once!
                    elementsOfShape = disassembleElsetToEnsightShapes(elSet)

                    thedata = evalModelAccessibleExpression(definition["modeldata"], modelInfo)
                    result = np.asarray(thedata, dtype=float)

                    varDict = {shape: result[elIndicesOfShape] for shape, elIndicesOfShape in elementsOfShape.items()}

                    if len(result.shape) == 1:
                        dimension = 1
                    else:
                        dimension = result.shape[1]

                    partsDict = {part.partNumber: varDict}
                    enSightVar = EnsightPerElementVariable(name, dimension, partsDict)
                    self.ensightCase.writeVariableTrendChunk(enSightVar, self.staticTAndFSetNumber)
                    del enSightVar

    def initializeStep(self, step, stepActions, stepOptions):
        if self.name in stepOptions or "Ensight" in stepOptions:
            options = stepOptions.get(self.name, False) or stepOptions["Ensight"]
            self.intermediateSaveInterval = int(options.get("intermediateSaveInterval", self.intermediateSaveInterval))
            self.minDTForOutput = float(options.get("minDTForOutput", self.minDTForOutput))

    def finalizeIncrement(self, U, P, increment):
        incNumber, incrementSize, stepProgress, dT, stepTimeAtIncrementStart, totalTimeAtIncrementStart = increment
        time = totalTimeAtIncrementStart + dT

        # check if we should write output, i.e., if enough time has passed:
        timeSinceLastOutput = time - self.timeAtLastOutput
        if  self.minDTForOutput - timeSinceLastOutput > 1e-12 :
            self.journal.message( "skipping output", self.identification, 1)
            return 

        if dT <= 0.0 and self.finishedSteps > 0:
            self.journal.message(
                "skipping output for zero-increment in step {:}".format(self.finishedSteps + 1), self.identification, 1
            )
            return

        self.writeOutput ( U, P, time)

    def writeOutput(self, U, P, time):

        self.timeAtLastOutput = time 
        self.ensightCase.setCurrentTime(self.transientTAndFSetNumber, time)

        for perNodeJob in self.transientPerNodeJobs:
            resultTypeLength = perNodeJob["varSize"]
            jobName = perNodeJob["name"]
            nodalVarTable = perNodeJob["fieldOutput"].getLastResult()
            partsDict = {perNodeJob["part"].partNumber: ("coordinates", nodalVarTable)}
            enSightVar = EnsightPerNodeVariable(jobName, resultTypeLength, partsDict)
            self.ensightCase.writeVariableTrendChunk(enSightVar, self.transientTAndFSetNumber)
            del enSightVar

        for perElementJob in self.transientPerElementJobs:
            name = perElementJob["name"]
            part = perElementJob["part"]
            elementsOfShape = perElementJob["elementsOfShape"]
            result = perElementJob["fieldOutput"].getLastResult()
            varDict = {shape: result[elIndicesOfShape] for shape, elIndicesOfShape in elementsOfShape.items()}

            if len(result.shape) == 1:
                dimension = 1
            else:
                dimension = result.shape[1]

            partsDict = {part.partNumber: varDict}
            enSightVar = EnsightPerElementVariable(name, dimension, partsDict)
            self.ensightCase.writeVariableTrendChunk(enSightVar, self.transientTAndFSetNumber)
            del enSightVar

        # intermediate save of the case
        if self.intermediateSaveInterval:
            if self.intermediateSaveIntervalCounter >= self.intermediateSaveInterval:
                self.ensightCase.finalize(replaceTimeValuesByEnumeration=False, closeFileHandes=False)
                self.intermediateSaveIntervalCounter = 0
            self.intermediateSaveIntervalCounter += 1

    def finalizeStep(
        self,
        U,
        P,
        time,
    ):
        if time - self.timeAtLastOutput > 1e-12:
            self.writeOutput(U, P, time)

        self.finishedSteps += 1

    def finalizeJob(
        self,
        U,
        P,
    ):
        self.ensightCase.finalize(replaceTimeValuesByEnumeration=False)
