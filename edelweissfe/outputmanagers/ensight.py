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
# Created on Sun Jan 15 14:22:48 2017

# @author: Matthias Neuner

import datetime
import os
from collections import defaultdict
from distutils.util import strtobool
from io import TextIOBase

import numpy as np

from edelweissfe.models.femodel import FEModel
from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from edelweissfe.points.node import Node
from edelweissfe.sets.elementset import ElementSet
from edelweissfe.sets.nodeset import NodeSet
from edelweissfe.utils.fieldoutput import (
    ElementFieldOutput,
    NodeFieldOutput,
    _FieldOutputBase,
)
from edelweissfe.utils.meshtools import disassembleElsetToEnsightShapes

"""
Output manager for Ensight exports.
If loaded, it automatically exports all elSets as Ensight parts.
For each part, perNode and perElement results can be exported, which are imported from fieldOutputs.

"""

documentation = {
    "configuration": "overwrite=y|n|yes|no|True|False|",
    "create": "The of the variabl; valid values: perNode|perElement",
    "fieldOutput": "Name of the result, defined on an elSet (also for perNode results!)",
    "name": "(Optional), default = the field output's name",
    "intermediateSaveInterval": 'Step option in category "Ensight": save .case file every N increments',
}


def writeCFloat(f, ndarray):
    np.asarray(ndarray, dtype=np.float32).tofile(f)


def writeCInt(f, ndarray):
    np.asarray(ndarray, dtype=np.int32).tofile(f)


def writeC80(f, string):
    np.asarray(string, dtype="S80").tofile(f)


ensightPerNodeVariableTypes = {
    1: "scalar per node",
    3: "vector per node",
    6: "tensor symm per node",
    9: "tensor asym per node",
}

ensightPerElementVariableTypes = {
    1: "scalar per element",
    3: "vector per element",
    6: "tensor symm per element",
    9: "tensor asym per element",
}


class EnsightUnstructuredPart:
    """Represents an unstructured ENsight part, defined by a list of nodes and a dictionary of elements.
    Each dictionary entry consists of a list of tuples of elementlabel and nodelist:

    Parameters
    ----------
    description
        A string describing the name of this part
    partNumber
        A unique integer identifying this part
    nodes
        The list of nodes in this part
    elementTree
        A dictionary, with entries for
            - each element shape,
            - containing a list of elements
                - defined by a tuple of
                    - a label and
                    - the node index list.
    """

    def __init__(
        self,
        description: str,
        partNumber: int,
        nodes: list[Node],
        elementTree: dict[str, dict[int, list[Node]]],
    ):
        self.structureType = "coordinates"
        self.elementTree = elementTree

        self.nodes = nodes
        self.nodeLabels = np.asarray([node.label for node in nodes], np.int32)

        self.description = description  # string, describing the part; max. 80 characters
        self.partNumber = partNumber
        self.nodeCoordinateArray = np.asarray([node.coordinates for node in nodes])

    def writeToFile(
        self,
        binaryFileHandle=TextIOBase,
        printNodeLabels: bool = True,
        printElementLabels: bool = True,
    ):
        """
        Write the part to a file.

        Parameters
        ----------
        binaryFileHandle
            The file handle for writing.
        printNodeLabels
            Write the node labels.
        printElementLabels
            Write the element labels.
        """

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
            writeCFloat(
                f,
                np.zeros(self.nodeCoordinateArray.shape[0] * (3 - self.nodeCoordinateArray.shape[1])),
            )

        for elemType, elements in self.elementTree.items():
            writeC80(f, elemType)
            writeCInt(f, len(elements))
            if printElementLabels:
                writeCInt(f, np.asarray([elNumber for elNumber in elements.keys()], np.int32))

            for nodeIndices in elements.values():
                writeCInt(f, np.asarray(nodeIndices, np.int32)[:] + 1)


class EnsightTimeSet:
    """Represents a set which may be used by EnsightGeometry, EnsightStructuredPart, EnsightUnstructuredPart and is written into the case file.

    Parameters
    ----------
    number
        The unique number of this time set.
    description
        A description of this time set.
    fileNameStartNumber
        Counter start of the multiple files counter.
    fileNameNumberIncrement
        Counter increment of the multpple files counter.
    timeValues
        A list of time values.
    """

    def __init__(
        self,
        number: int = 1,
        description: str = "timeStepDesc",
        fileNameStartNumber: int = 0,
        fileNameNumberIncrement: int = 1,
        timeValues: list = None,
    ):
        self.number = number
        self.description = description
        self.fileNameStartNumber = fileNameStartNumber
        self.fileNameNumberIncrement = fileNameNumberIncrement
        self.timeValues = timeValues if timeValues is not None else []


class EnsightGeometry:
    """Container class for one or more EnsightParts at a certain time state, handles also the file writing operation.

    Parameters
    ----------
    name
        The name.
    descriptionLine1
        Description line 1.
    descriptionLine2
        Description line 2.
    ensightPartList
        The list of unstructured parts in this geometry.
    nodeIdOption
        Ensight option for handling the node ids.
    elementIdOption
        Ensight option for handling the elemnt ids.
    """

    def __init__(
        self,
        name: str = "geometry",
        descriptionLine1: str = "",
        descriptionLine2: str = "",
        ensightPartList: list[EnsightUnstructuredPart] = None,
        nodeIdOption: str = "given",
        elementIdOption: str = "given",
    ):
        self.name = name
        self.descLine1 = descriptionLine1
        self.descLine2 = descriptionLine2
        self.partList = ensightPartList if ensightPartList is not None else []
        self.nodeIdOption = nodeIdOption
        self.elementIdOption = elementIdOption

    def writeToFile(self, fileHandle: TextIOBase):
        """Write the variable to a file."

        Parameters
        ----------
        fileHandle
            The file handle for writing the file.
        """
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


class EnsightGeometryTrend:
    """Container class for the time dependent evolution of the geometry,
    establishes the connection between the geometry entities and a EnsightTimeSet.

    Parameters
    ----------
    ensightTimeSet
        The Timeset.
    ensightGeometryList
        A list of evolving Ensight Geometries.
    """

    def __init__(
        self,
        ensightTimeSet: EnsightTimeSet,
        ensightGeometryList: list[EnsightGeometry] = None,
    ):
        self.timeSet = ensightTimeSet
        self.geometryList = ensightGeometryList if ensightGeometryList is not None else []


class EnsightVariableTrend:
    """Container class for the time dependent evolution of one variable,
    establishes the connection between EnsightVariable entities and a EnsighTimeSet.

    Parameters
    ----------
    ensightTimeSet
        The timeset.
    variableName
        The name of this variable.
    ensightVariableList
        The list of variables over time.
    variableType
        The Ensight valid type of this variable.
    description
        The description of this variable.
    """

    def __init__(
        self,
        ensightTimeSet: EnsightTimeSet,
        variableName: str,
        ensightVariableList: list = None,
        variableType="scalar per node",
        description="variableTrendDescription",
    ):
        self.timeSet = ensightTimeSet
        self.variableName = variableName
        self.variableList = ensightVariableList if ensightVariableList is not None else []
        self.variableType = variableType
        self.description = description


class EnsightPerNodeVariable:
    """Container class for data for one certain variable, defined for one or more parts (classification by partID), at a certain time state.
    For each part the structuretype ("coordinate" or "block") has to be defined.
    Each part-variable assignment is defined by a dictionary entry of type: { EnsightPart: np.array(variableValues) }

    Parameters
    ----------
    name
        The name of this variable.
    variableDimension
        The size of the variable per value.
    ensightPartsDict
        A dictionary defining the values for given Ensight parts.
    """

    def __init__(
        self,
        name: str,
        ensightPartsDict: dict[EnsightUnstructuredPart, np.ndarray],
        varSize: int,
    ):

        self.name = name.replace(" ", "_")
        self.description = self.name
        self.partsDict = ensightPartsDict or {}  # { EnsightPart: np.array(variableValues) }

        self.varType = ensightPerNodeVariableTypes[varSize]

    def writeToFile(
        self,
        fileHandle: TextIOBase,
    ):
        """Write the variable to a file.

        Parameters
        ----------
        fileHandle
            The file handle for writing.
        """

        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, (structureType, values) in self.partsDict.items():
            writeC80(f, "part")
            writeCInt(f, ensightPartID)
            writeC80(f, structureType)
            writeCFloat(f, values.T)


class EnsightPerElementVariable:
    """Container class for data for one certain variable, defined for one or more parts (classification by partID), at a certain time state.
    For each part the structuretype ("coordinate" or "block") has to be defined.
    Each part-variable assignment is defined by a dictionary entry of type: { EnsightPart: np.array(variableValues) }

    Parameters
    ----------
    name
        The name of this variable.
    variableDimension
        The size of this variable per entry.
    ensightPartsDict
        The dictionary containing parts and their values.
    """

    def __init__(
        self,
        name: str,
        ensightPartsDict: dict[EnsightUnstructuredPart, np.ndarray],
        varSize: int,
    ):
        self.name = name.replace(" ", "_")
        self.description = self.name
        self.partsDict = ensightPartsDict
        self.varType = ensightPerElementVariableTypes[varSize]

    def writeToFile(self, fileHandle: TextIOBase):
        """Write the variable to a file.

        Parameters
        ----------
        fileHandle
            The file handle for writing.
        """

        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, elTypeDict in self.partsDict.items():
            writeC80(f, "part")
            writeCInt(f, ensightPartID)
            for elType, values in elTypeDict.items():
                writeC80(f, elType)
                writeCFloat(f, values.T)


class EnsightChunkWiseCase:
    """An Ensight case, containg of time sets, geometry trends and variable trends,
    which can be written in chunks at certain times.

    Parameters
    ----------
    caseName
        The name of this case.
    directory
        The path to write.
    writeTransientSingleFiles
        Write single or multiple files.
    """

    def __init__(self, caseName: str, directory: str = "", writeTransientSingleFiles: bool = True):
        self.directory = directory
        self.caseName = caseName
        self.caseFileNamePrefix = os.path.join(directory, caseName)
        self.writeTransientSingleFiles = writeTransientSingleFiles
        self.timeAndFileSets = {}
        self.geometryTrends = {}
        self.variableTrends = {}
        self.fileNames = {}

        if not os.path.exists(self.caseFileNamePrefix):
            os.mkdir(self.caseFileNamePrefix)

    def setCurrentTime(self, timeAndFileSetNumber: int, timeValue: float):
        """Set the current time of the case.
        Parameters
        ----------
        timeAndFileSetNumber
            The number of the file and time set.
        timeValue
            The time value.
        """

        if timeAndFileSetNumber not in self.timeAndFileSets:
            self.timeAndFileSets[timeAndFileSetNumber] = EnsightTimeSet(timeAndFileSetNumber, "no description", 0, 1)
        tfSet = self.timeAndFileSets[timeAndFileSetNumber]
        tfSet.timeValues.append(timeValue)

    def writeGeometryTrendChunk(self, ensightGeometry: EnsightGeometryTrend, timeAndFileSetNumber: int = 1):
        """
        Write a chunk of geometry trend.

        Parameters
        ----------
        ensightGeometry
            The trend to write.
        timeAndFileSetNumber
            The associated time and fileset number.
        """

        if ensightGeometry.name not in self.fileNames:
            fileName = os.path.join(
                self.caseFileNamePrefix,
                ensightGeometry.name + ".geo",
            )

            self.fileNames[ensightGeometry.name] = fileName
            # create empty file
            with open(fileName, mode="wb") as f:
                pass

        filename = self.fileNames[ensightGeometry.name]

        with open(filename, mode="ab") as f:
            if ensightGeometry.name not in self.geometryTrends:
                self.geometryTrends[ensightGeometry.name] = timeAndFileSetNumber
                writeC80(f, "C Binary")

            if self.writeTransientSingleFiles:
                writeC80(f, "BEGIN TIME STEP")
                ensightGeometry.writeToFile(f)
                writeC80(f, "END TIME STEP")

    def writeVariableTrendChunk(self, ensightVariable: EnsightVariableTrend, timeAndFileSetNumber: int = 2):
        """
        Write a chunk of variable trend.

        Parameters
        ----------
        ensightVariable
            The trend to write.
        timeAndFileSetNumber
            The associated time and fileset number.
        """

        if ensightVariable.name not in self.fileNames:
            # create file name
            fileName = os.path.join(self.caseFileNamePrefix, ensightVariable.name + ".var")

            # append to file names
            self.fileNames[ensightVariable.name] = fileName

            # create empty file
            with open(fileName, mode="wb") as f:
                pass

        filename = self.fileNames[ensightVariable.name]

        with open(filename, mode="ab") as f:

            if ensightVariable.name not in self.variableTrends:
                self.variableTrends[ensightVariable.name] = (
                    timeAndFileSetNumber,
                    ensightVariable.varType,
                )
                writeC80(f, "C Binary")

            if self.writeTransientSingleFiles:
                writeC80(f, "BEGIN TIME STEP")
                ensightVariable.writeToFile(f)
                writeC80(f, "END TIME STEP")

    def finalize(self, replaceTimeValuesByEnumeration: bool = True, closeFileHandes: bool = True):
        """Write the file .case file containing all the required information."

        Parameters
        ----------
        replaceTimeValuesByEnumeration
            Remove the factor of time and make integers as discrete steps only.
        closeFileHandes
            Close the file handles after writing.
        """

        caseFName = self.caseFileNamePrefix + ".case"

        with open(caseFName, mode="w") as cf:
            cf.write("FORMAT\n")
            cf.write("type: ensight gold\n")

            cf.write("TIME\n")
            for setNum, timeSet in self.timeAndFileSets.items():
                cf.write("time set: " + str(setNum) + " no description\n")
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
                    cf.write("file set: {:}\n".format(timeSet.number))
                    cf.write("number of steps: {:}\n".format(len(timeSet.timeValues)))

            cf.write("GEOMETRY\n")
            for geometryName, tAndFSetNum in self.geometryTrends.items():
                cf.write(
                    "model: {:} {:} {:}\n".format(
                        tAndFSetNum,
                        tAndFSetNum,
                        os.path.join(self.caseFileNamePrefix, geometryName + ".geo"),
                    )
                )

            cf.write("VARIABLE\n")
            for variableName, (
                tAndFSetNum,
                variableType,
            ) in self.variableTrends.items():
                cf.write(
                    "{:}: {:} {:} {:} {:}.var\n".format(
                        variableType,
                        tAndFSetNum,
                        tAndFSetNum,
                        variableName,
                        os.path.join(self.caseFileNamePrefix, variableName),
                    )
                )


def createUnstructuredPartFromElementSet(setName, elementSet: list, partID: int):
    """Determines the element and node list for an Ensightpart from an
    element set. The reduced, unique node set is generated, as well as
    the element to node index mapping for the ensight part.

    Parameters
    ----------
    elementSet
        The list of elements defining this part.
    partID
        The id of this part.
    """

    nodeCounter = 0
    partNodes = dict()
    elementDict = dict()
    for element in elementSet:
        elShape = element.ensightType
        if elShape not in elementDict:
            elementDict[elShape] = dict()
        elNodeIndices = []
        for node in element.visualizationNodes:
            # if the node is already in the dict, get its index,
            # else insert it, and get the current idx = counter. increase the counter
            idx = partNodes.setdefault(node, nodeCounter)
            elNodeIndices.append(idx)
            if idx == nodeCounter:
                # the node was just inserted, so increase the counter of inserted nodes
                nodeCounter += 1
        elementDict[elShape][element.elNumber] = elNodeIndices

    return EnsightUnstructuredPart(setName, partID, partNodes.keys(), elementDict)


def createUnstructuredPartFromNodeSet(setName, nodeSet: list, partID: int):
    """Construtcts an EnSight part for a node set. Since EnSight parts comprise nodes and elements, each node is assigned to a dummy point element.

    Parameters
    ----------
    nodeSet
        The list of nodes defining this part.
    partID
        The id of this part.
    """

    elementDict = dict()

    elementDict["point"] = {i: [i] for i in range(len(nodeSet))}

    return EnsightUnstructuredPart("NSET_" + setName, partID, list(nodeSet), elementDict)


class OutputManager(OutputManagerBase):
    identification = "Ensight Export"

    def __init__(self, name, model, fieldOutputController, journal, plotter):
        self.name = name

        self.model = model
        self.timeAtLastOutput = -1e16
        self.minDTForOutput = -1e16
        self.finishedSteps = 0
        self.intermediateSaveInterval = 10
        self.intermediateSaveIntervalCounter = 0
        self.fieldOutputController = fieldOutputController
        self.journal = journal

        self.transientTAndFSetNumber = 1
        self.staticTAndFSetNumber = 2

        self.elSetToEnsightPartMappings = {}
        self.nSetToEnsightPartMappings = {}

        self._transientPerNodeVariableJobs = defaultdict(list)
        self._transientPerElementVariableJobs = defaultdict(list)

        self.exportName = name

        self.geometryParts = self._createGeometryParts(1)

    def updateDefinition(self, **kwargs: dict):
        self.model
        # standard, transient jobs accessing the fieldoutput:

        # Determine the type
        if "create" in kwargs:
            create = kwargs.pop("create")
            fieldOutput = kwargs.pop("fieldOutput")
            part = None
            if "nSet" in kwargs:
                part = self.nSetToEnsightPartMappings[kwargs.pop("nSet")]
            elif "elSet" in kwargs:
                part = self.elSetToEnsightPartMappings[kwargs.pop("elSet")]

            name = kwargs.get("name", fieldOutput.name).replace(" ", "_")

            nEntries, varSize = self._ensureArrayIs2D(fieldOutput.getLastResult()).shape

            if self.model.domainSize == 2 and varSize == 2:
                varSize = 3

            transient = kwargs.get("transient", "True")
            transient = strtobool(transient)

            if create == "perElement":
                self.createPerElementOutput(fieldOutput, part, name, transient=transient, varSize=varSize)
            elif create == "perNode":
                self.createPerNodeOutput(fieldOutput, part, name, transient=transient, varSize=varSize)

        if "configuration" in kwargs:
            # ensight output is overwritten by default
            self.overwrite = strtobool(kwargs.get("overwrite", "True"))
            if not self.overwrite:
                self.exportName = "{:}_{:}".format(self.name, datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))

    def createPerElementOutput(
        self,
        fieldOutput: ElementFieldOutput,
        part: set = None,
        name: str = None,
        transient: bool = True,
        varSize: int = None,
    ):
        """Create a per element output job.

        Parameters
        ----------
        fieldOutput
            The field output to export.
        part
            The part to which the output belongs. If not specified, the part is determined from the field output.
        name
            The name of the output. If not specified, the name of the field output is taken.
        transient
            Whether the output is transient.
        varSize
            The size of the variable. If not specified, the size of the field output is taken.
        """

        variableJob = dict()
        variableJob["fieldOutput"] = fieldOutput
        variableJob["part"] = part
        if not name:
            name = fieldOutput.name
        variableJob["name"] = name
        variableJob["transient"] = transient

        nEntries, varSizeFp = self._ensureArrayIs2D(fieldOutput.getLastResult()).shape

        if not part:
            part = self._getTargetPartForFieldOutput(fieldOutput)
        variableJob["part"] = part

        if nEntries != len(fieldOutput.associatedSet):
            raise Exception(
                "Variable {:} result size ({:}) does not match the number of nodes ({:})".format(
                    variableJob["name"], nEntries, len(variableJob["part"].nodes)
                )
            )

        if not varSize:
            varSize = varSizeFp
        variableJob["varSize"] = varSize

        variableJob["elementsOfShape"] = disassembleElsetToEnsightShapes(fieldOutput.associatedSet)

        if transient:
            self._transientPerElementVariableJobs[variableJob["name"]].append(variableJob)
        else:
            raise Exception("Only transient per element outputs are supported!")

    def createPerNodeOutput(
        self,
        fieldOutput: NodeFieldOutput,
        part: EnsightUnstructuredPart = None,
        name: str = None,
        transient: bool = True,
        varSize: int = None,
    ):
        """Create a per node output job.

        Parameters
        ----------
        fieldOutput
            The field output to export.
        part
            The part to which the output belongs. If not specified, the part is determined from the field output.
        name
            The name of the output. If not specified, the name of the field output is taken.
        transient
            Whether the output is transient.
        varSize
            The size of the variable. If not specified, the size of the field output is taken.
        """

        variableJob = dict()
        variableJob["fieldOutput"] = fieldOutput
        variableJob["part"] = part
        if not name:
            name = fieldOutput.name
        variableJob["name"] = name
        variableJob["transient"] = transient

        nEntries, varSizeFp = self._ensureArrayIs2D(fieldOutput.getLastResult()).shape
        if not varSize:
            varSize = varSizeFp

        variableJob["varSize"] = varSize

        if not part:
            part = self._getTargetPartForFieldOutput(fieldOutput)
        variableJob["part"] = part

        if nEntries != len(variableJob["part"].nodes):
            raise Exception(
                "Variable {:} result size ({:}) does not match the number of nodes ({:})".format(
                    variableJob["name"], nEntries, len(variableJob["part"].nodes)
                )
            )

        if transient:
            self._transientPerNodeVariableJobs[variableJob["name"]].append(variableJob)
        else:
            raise Exception("Only transient per node outputs are supported!")

    def initializeJob(self):
        self.ensightCase = EnsightChunkWiseCase(self.exportName)
        self.ensightCase.setCurrentTime(self.staticTAndFSetNumber, self.model.time)

        geometry = EnsightGeometry("geometry", "EdelweissFE", "*export*", ensightPartList=self.geometryParts)
        geometryTimesetNumber = self.staticTAndFSetNumber

        self.ensightCase.writeGeometryTrendChunk(geometry, geometryTimesetNumber)

    def initializeStep(self, step):
        if self.name in step.actions["options"] or "Ensight" in step.actions["options"]:
            options = step.actions["options"].get(self.name, False) or step.actions["options"]["Ensight"].options
            self.intermediateSaveInterval = int(options.get("intermediateSaveInterval", self.intermediateSaveInterval))
            self.minDTForOutput = float(options.get("minDTForOutput", self.minDTForOutput))

    def finalizeIncrement(self, **kwargs):
        time = self.model.time

        # check if we should write output, i.e., if enough time has passed:
        timeSinceLastOutput = time - self.timeAtLastOutput

        if self.minDTForOutput - timeSinceLastOutput > 1e-16:
            self.journal.message("Skipping output".format(), self.identification, 1)
            return

        self.writeOutput(self.model)

    def finalizeFailedIncrement(self, **kwargs):
        pass

    def writeOutput(self, model: FEModel):
        self.timeAtLastOutput = model.time
        self.ensightCase.setCurrentTime(self.transientTAndFSetNumber, model.time)

        resultsByParts = {}
        for resultName, perNodeVariableJobs in self._transientPerNodeVariableJobs.items():
            for perNodeVariableJob in perNodeVariableJobs:
                result = self._ensureArrayIs2D(perNodeVariableJob["fieldOutput"].getLastResult())

                if self.model.domainSize == 2 and result.shape[1] == 2:
                    result = self._make2DVector3D(result)

                resultsByParts[perNodeVariableJob["part"].partNumber] = ("coordinates", result)
            enSightVariable = EnsightPerNodeVariable(resultName, resultsByParts, perNodeVariableJob["varSize"])
            self.ensightCase.writeVariableTrendChunk(enSightVariable, self.transientTAndFSetNumber)
            del enSightVariable

        for resultName, perElementVariableJobs in self._transientPerElementVariableJobs.items():
            resultsByParts = {}
            for perElementVariableJob in perElementVariableJobs:
                part = perElementVariableJob["part"]
                elementsOfShape = perElementVariableJob["elementsOfShape"]
                result = self._ensureArrayIs2D(perElementVariableJob["fieldOutput"].getLastResult())

                if self.model.domainSize == 2 and result.shape[1] == 2:
                    result = self._make2DVector3D(result)

                partResultsByElementShape = {
                    shape: result[elIndicesOfShape] for shape, elIndicesOfShape in elementsOfShape.items()
                }
                resultsByParts[part.partNumber] = partResultsByElementShape

            enSightVariable = EnsightPerElementVariable(resultName, resultsByParts, perElementVariableJob["varSize"])
            self.ensightCase.writeVariableTrendChunk(enSightVariable, self.transientTAndFSetNumber)
            del enSightVariable

        # intermediate save of the case
        if self.intermediateSaveInterval:
            if self.intermediateSaveIntervalCounter >= self.intermediateSaveInterval:
                self.ensightCase.finalize(replaceTimeValuesByEnumeration=False, closeFileHandes=False)
                self.intermediateSaveIntervalCounter = 0
            self.intermediateSaveIntervalCounter += 1

    def finalizeStep(
        self,
    ):
        if self.model.time - self.timeAtLastOutput > 1e-12:
            self.writeOutput(self.model)

        self.finishedSteps += 1

    def finalizeJob(
        self,
    ):
        self.ensightCase.finalize(replaceTimeValuesByEnumeration=False)

    def _createGeometryParts(self, firstPartID: int):
        model = self.model
        elementSets = model.elementSets

        elSetParts = []
        partCounter = firstPartID
        for setName, elSet in elementSets.items():
            elSetPart = createUnstructuredPartFromElementSet(setName, elSet, partCounter)
            self.elSetToEnsightPartMappings[setName] = elSetPart
            elSetParts.append(elSetPart)
            partCounter += 1

        nodeSets = model.nodeSets

        nodeSetParts = []
        for setName, nodeSet in nodeSets.items():
            nodeSetPart = createUnstructuredPartFromNodeSet(setName, nodeSet, partCounter)
            self.nSetToEnsightPartMappings[setName] = nodeSetPart
            nodeSetParts.append(nodeSetPart)
            partCounter += 1

        return elSetParts + nodeSetParts

    def _getTargetPartForFieldOutput(self, fieldOutput: _FieldOutputBase) -> EnsightUnstructuredPart:
        """
        Determine depending on the specified input,
        for which Part the result should be written.
        If no input is specified, the associated set of the :class:`FieldOutput is taken.

        Parameters
        ----------
        fieldOutput
            The :class:`FieldOutput which contains the result

        Returns
        -------
        EnsightStructuredPart
            The identified part.
        """

        theSetName = fieldOutput.associatedSet.name

        if isinstance(fieldOutput.associatedSet, NodeSet):
            return self.nSetToEnsightPartMappings[theSetName]

        elif isinstance(fieldOutput.associatedSet, ElementSet):
            return self.elSetToEnsightPartMappings[theSetName]

        else:
            raise Exception(
                "Ensight Variables need to be excplicity associated with a part, our implicitly through a FieldOutput defined on ElementSets or NodeSets!"
            )

    def _ensureArrayIs2D(self, result: np.ndarray) -> np.ndarray:
        """
        Ensure that a result array is in tabular form.

        Parameters
        ----------
        result
            The result in potential 1D form.

        Returns
        -------
        np.ndarray
            The result in guaranteed 2d form.
        """
        return np.reshape(result, (len(result), -1))

    def _make2DVector3D(self, result: np.ndarray) -> np.ndarray:
        """
        Vector results in 2D consist of 2 components.
        However, in Ensight all results must be written for the 3D case.
        This function appends a zero column.

        Parameters
        ----------
        result
            The vector result in 2D form.

        Returns
        -------
        np.ndarray
            The result in 3D form.
        """

        return np.pad(
            result,
            ((0, 0), (0, 1)),
        )
