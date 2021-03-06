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
Created on Tue Apr  3 09:44:10 2018

@author: Matthias Neuner
"""

from fe.elements.node import Node
from fe.config.elementlibrary import getElementByName
from fe.utils.misc import isInteger
from fe.config.constraints import getConstraintByName
from fe.config.sections import getSectionByName
from fe.config.analyticalFields import getAnalyticalFieldByName
import numpy as np


class AbqModelConstructor:
    def __init__(self, journal):
        self.journal = journal

    def createGeometryFromInputFile(self, modelInfo, inputFile):
        """Collects nodes, elements, node sets and element sets from
        the input file."""

        domainSize = modelInfo["domainSize"]

        # returns an dict of {node label: node}
        nodeDefinitions = modelInfo["nodes"]
        for nodeDefs in inputFile["*node"]:
            currNodeDefs = {}
            for defLine in nodeDefs["data"]:
                label = int(defLine[0])
                coordinates = np.zeros(domainSize)
                coordinates[:] = defLine[1:]
                currNodeDefs[label] = Node(
                    label,
                    coordinates,
                )
            nodeDefinitions.update(currNodeDefs)

            if "nset" in nodeDefs.keys():
                setName = nodeDefs["nset"]
                modelInfo["nodeSets"][setName] = list(currNodeDefs.keys())

        # returns an dict of {element Label: element}
        elements = modelInfo["elements"]

        for elDefs in inputFile["*element"]:
            elementType = elDefs["type"]
            elementProvider = elDefs.get("provider", False)
            ElementClass = getElementByName(elementType, elementProvider)

            currElDefs = {}
            for defLine in elDefs["data"]:
                label = defLine[0]
                # store nodeObjects in elNodes list
                elNodes = [nodeDefinitions[n] for n in defLine[1:]]
                newEl = ElementClass(elementType, label)
                newEl.setNodes(elNodes)
                currElDefs[label] = newEl
            elements.update(currElDefs)

            if "elset" in elDefs.keys():
                setName = elDefs["elset"]
                modelInfo["elementSets"][setName] = list(currElDefs.values())

        # generate dictionary of elementObjects belonging to a specified elementset
        # or generate elementset by generate definition in inputfile
        elementSets = modelInfo["elementSets"]

        for elSetDefinition in inputFile["*elSet"]:
            name = elSetDefinition["elSet"]

            # decide if entries are labels or existing nodeSets:
            if isInteger(elSetDefinition["data"][0][0]):
                elNumbers = [int(num) for line in elSetDefinition["data"] for num in line]

                if elSetDefinition.get("generate", False):
                    generateDef = elNumbers[0:3]
                    els = [
                        elements[n] for n in np.arange(generateDef[0], generateDef[1] + 1, generateDef[2], dtype=int)
                    ]

                elif elSetDefinition.get("boolean", False):
                    booleanDef = elSetDefinition.get("boolean")
                    if booleanDef == "difference":
                        els = [n for n in elementSets[name] if n.elNumber not in elNumbers]

                    elif booleanDef == "union":
                        els = [n for n in elementSets[name]]
                        els += [elements[n] for n in elNumbers]

                    elif booleanDef == "intersection":
                        elNumbersBase = [n.elNumber for n in elementSets[name]]
                        els = [elements[n] for n in list(set(elNumbers).intersection(elNumbersBase))]
                    else:
                        raise Exception("Undefined boolean operation!")

                    if elSetDefinition.get("newElSet") != name:
                        name = elSetDefinition.get("newElSet")
                    else:
                        del elementSets[name]
                else:
                    els = [elements[elNum] for elNum in elNumbers]
                elementSets[name] = els
            else:
                elementSets[name] = []
                for line in elSetDefinition["data"]:
                    for elSet in line:
                        elementSets[name] += elementSets[elSet]

        # generate dictionary of nodeObjects belonging to a specified nodeset
        # or generate nodeset by generate definition in inputfile
        nodeSets = modelInfo["nodeSets"]
        for nSetDefinition in inputFile["*nSet"]:
            name = nSetDefinition["nSet"]

            if isInteger(nSetDefinition["data"][0][0]):
                nodes = [int(n) for line in nSetDefinition["data"] for n in line]
                if nSetDefinition.get("generate", False):
                    generateDef = nodes  # nSetDefinition['data'][0][0:3]
                    nodes = [
                        nodeDefinitions[n]
                        for n in np.arange(generateDef[0], generateDef[1] + 1, generateDef[2], dtype=int)
                    ]
                else:
                    nodes = [nodeDefinitions[n] for n in nodes]
                nodeSets[name] = nodes
            else:
                nodeSets[name] = []
                for line in nSetDefinition["data"]:
                    for nSet in line:
                        nodeSets[name] += nodeSets[nSet]

        modelInfo["nodeSets"]["all"] = list(modelInfo["nodes"].values())
        modelInfo["elementSets"]["all"] = list(modelInfo["elements"].values())

        # generate surfaces sets
        for surfaceDef in inputFile["*surface"]:
            name = surfaceDef["name"]
            sType = surfaceDef.get("type", "element").lower()
            surface = {}
            if sType == "element":
                for l in surfaceDef["data"]:
                    elSet, faceNumber = l
                    faceNumber = int(faceNumber.replace("S", ""))
                    elements = modelInfo["elementSets"][elSet]
                    elements += surface.setdefault(faceNumber, [])
                    surface[faceNumber] = elements

            modelInfo["surfaces"][name] = surface

        return modelInfo

    def createMaterialsFromInputFile(self, modelInfo, inputFile):
        for materialDef in inputFile["*material"]:

            materialName = materialDef["name"]
            materialID = materialDef.get("id", materialName)
            materialProperties = np.hstack(materialDef["data"])

            modelInfo["materials"][materialID] = {"name": materialName, "properties": materialProperties}

        return modelInfo

    def createConstraintsFromInputFile(self, modelInfo, inputFile):
        for constraintDef in inputFile["*constraint"]:
            name = constraintDef["name"]
            constraint = constraintDef["type"]
            data = constraintDef["data"]

            constraint = getConstraintByName(constraint)(name, data, modelInfo)
            modelInfo["constraints"][name] = constraint

        return modelInfo

    def createSectionsFromInputFile(self, modelInfo, inputFile):
        """Assign properties and section properties to all elements by
        the given section definitions."""

        elementSets = modelInfo["elementSets"]

        for secDef in inputFile["*section"]:
            name = secDef["name"]
            sec = secDef["type"]
            data = secDef["data"]
            materialID = secDef["material"]

            Section = getSectionByName(sec)

            # this was a bad design decision, and will be deprecated sooner or later:
            thickness = float(secDef.get("thickness", 0.0))

            theSection = Section(name, data, materialID, thickness, modelInfo)

            modelInfo = theSection.assignSectionPropertiesToModel(modelInfo)

        return modelInfo

    def createAnalyticalFieldsFromInputFile(self, modelInfo, inputFile):
        for fieldDef in inputFile["*analyticalField"]:
            analyticalFieldName = fieldDef["name"]
            analyticalFieldType = fieldDef["type"]
            analyticalFieldData = fieldDef["data"]

            analyticalFieldClass = getAnalyticalFieldByName(analyticalFieldType)
            analyticalField = analyticalFieldClass(analyticalFieldName, analyticalFieldData, modelInfo)

            modelInfo["analyticalFields"][analyticalFieldName] = analyticalField

        return modelInfo
