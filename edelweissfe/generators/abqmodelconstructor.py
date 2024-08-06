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
# Created on Tue Apr  3 09:44:10 2018

# @author: Matthias Neuner
"""
The default way to create finite element meshes
is using the keywords

 * ``*node``
 * ``*element``
 * ``*nset``
 * ``*elset``
 * ``*surface``

employing an Abaqus-like syntax.
"""

import numpy as np

from edelweissfe.config.analyticalfields import getAnalyticalFieldByName
from edelweissfe.config.constraints import getConstraintClass
from edelweissfe.config.elementlibrary import getElementClass
from edelweissfe.config.materiallibrary import getMaterialClass
from edelweissfe.config.sections import getSectionClass
from edelweissfe.models.femodel import FEModel
from edelweissfe.points.node import Node
from edelweissfe.sets.elementset import ElementSet
from edelweissfe.sets.nodeset import NodeSet
from edelweissfe.utils.misc import convertLinesToFlatArray, isInteger, splitLineAtCommas


class AbqModelConstructor:
    def __init__(self, journal):
        pass

    def createGeometryFromInputFile(self, model: FEModel, inputFile: dict) -> dict:
        """Collects nodes, elements, node sets and element sets from
        the input file.

        Parameters
        ----------
        model
            A dictionary containing the model tree.
        inputFile
            A dictionary contaning the input file tree.

        Returns
        -------
        dict
            The updated model tree.
        """

        domainSize = model.domainSize

        # returns an dict of {node label: node}
        nodeDefinitions = model.nodes
        for nodeDefs in inputFile["*node"]:
            currNodeDefs = {}
            for line in nodeDefs["data"]:
                defLine = splitLineAtCommas(line)

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
                model.nodeSets[setName] = NodeSet(setName, [nodeDefinitions[x] for x in currNodeDefs.keys()])

        # returns an dict of {element Label: element}
        elements = model.elements

        for elDefs in inputFile["*element"]:
            elementType = elDefs["type"]
            elementProvider = elDefs.get("provider")
            ElementClass = getElementClass(elementProvider)

            currElDefs = {}
            for line in elDefs["data"]:
                defLine = [int(i) for i in splitLineAtCommas(line)]

                label = defLine[0]
                # store nodeObjects in elNodes list
                elNodes = [nodeDefinitions[n] for n in defLine[1:]]
                newEl = ElementClass(elementType, label)
                newEl.setNodes(elNodes)
                currElDefs[label] = newEl
            elements.update(currElDefs)

            if "elset" in elDefs.keys():
                setName = elDefs["elset"]
                model.elementSets[setName] = ElementSet(setName, currElDefs.values())

        # generate dictionary of elementObjects belonging to a specified elementset
        # or generate elementset by generate definition in inputfile
        elementSets = model.elementSets

        for elSetDefinition in inputFile["*elSet"]:
            name = elSetDefinition["elSet"]

            data = [splitLineAtCommas(line) for line in elSetDefinition["data"]]
            # decide if entries are labels or existing nodeSets:
            if isInteger(data[0][0]):
                elNumbers = [int(num) for line in data for num in line]

                if elSetDefinition.get("generate", False):
                    generateDef = elNumbers[0:3]
                    els = [
                        elements[n]
                        for n in np.arange(
                            generateDef[0],
                            generateDef[1] + 1,
                            generateDef[2],
                            dtype=int,
                        )
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
                elementSets[name] = ElementSet(name, set(els))
            else:
                elementSets[name] = []
                for line in data:
                    for elSet in line:
                        elementSets[name] += elementSets[elSet]

        # generate dictionary of nodeObjects belonging to a specified nodeset
        # or generate nodeset by generate definition in inputfile
        nodeSets = model.nodeSets
        for nSetDefinition in inputFile["*nSet"]:
            name = nSetDefinition["nSet"]

            data = [splitLineAtCommas(line) for line in nSetDefinition["data"]]
            if isInteger(data[0][0]):
                nodes = [int(n) for line in data for n in line]
                if nSetDefinition.get("generate", False):
                    generateDef = nodes  # nSetDefinition['data'][0][0:3]
                    nodes = [
                        nodeDefinitions[n]
                        for n in np.arange(
                            generateDef[0],
                            generateDef[1] + 1,
                            generateDef[2],
                            dtype=int,
                        )
                    ]
                else:
                    nodes = [nodeDefinitions[n] for n in nodes]
                nodeSets[name] = NodeSet(name, nodes)
            else:
                nSetLabels = [nSetLabel for line in data for nSetLabel in line]
                nodes = [n for nSetLabel in nSetLabels for n in nodeSets[nSetLabel].nodes]
                nodeSets[name] = NodeSet(name, nodes)

        model.nodeSets["all"] = NodeSet("all", model.nodes.values())
        model.elementSets["all"] = ElementSet("all", model.elements.values())

        # generate surfaces sets
        for surfaceDef in inputFile["*surface"]:
            name = surfaceDef["name"]
            sType = surfaceDef.get("type", "element").lower()
            surface = {}
            if sType == "element":
                data = [splitLineAtCommas(line) for line in surfaceDef["data"]]
                for line in data:
                    elSet, faceNumber = line
                    faceNumber = int(faceNumber.replace("S", ""))
                    surface[faceNumber] = model.elementSets[elSet]

            model.surfaces[name] = surface

        return model

    def createMaterialsFromInputFile(self, model, inputFile):
        """Collects material defintions from the input file.
        Creates instances of materials.

        Parameters
        ----------
        model
            A dictionary containing the model tree.
        inputFile
            A dictionary contaning the input file tree.

        Returns
        -------
        dict
            The updated model tree.
        """

        for materialDef in inputFile["*material"]:
            materialName = materialDef["name"]
            materialProvider = materialDef.get("provider", None)
            materialID = materialDef.get("id", materialName)

            materialProperties = convertLinesToFlatArray(materialDef["data"], dtype=float)
            materialClass = getMaterialClass(materialName, materialProvider)

            if materialClass is None:  # for Marmot
                model.materials[materialID] = {
                    "name": materialName,
                    "properties": materialProperties,
                }
            else:  # for DisplacementElement
                model.materials[materialID] = materialClass(materialProperties)

        return model

    def createConstraintsFromInputFile(self, model, inputFile):
        """Collects constraint defintions from the input file.

        Parameters
        ----------
        model
            A dictionary containing the model tree.
        inputFile
            A dictionary contaning the input file tree.

        Returns
        -------
        dict
            The updated model tree.
        """

        for constraintDef in inputFile["*constraint"]:
            name = constraintDef["name"]
            constraint = constraintDef["type"]
            data = constraintDef["data"]

            constraint = getConstraintClass(constraint)(name, data, model)
            model.constraints[name] = constraint

        return model

    def createSectionsFromInputFile(self, model, inputFile):
        """Collects section defintions from the input file.
        Assigns properties and section properties to all elements by
        the given section definitions.

        Parameters
        ----------
        model
            A dictionary containing the model tree.
        inputFile
            A dictionary contaning the input file tree.

        Returns
        -------
        dict
            The updated model tree.
        """

        for definition in inputFile["*section"]:
            try:
                name = definition.pop("name")
            except KeyError:
                raise KeyError(f"No name specified for section {name}.")
            try:
                sectionType = definition.pop("type")
            except KeyError:
                raise KeyError(f"No type specified for section {name}.")
            try:  # should data be required?
                data = definition.pop("data")
            except KeyError:
                raise KeyError(f"No data specified for section {name}.")
            try:
                materialName = definition.pop("material")
            except KeyError:
                raise KeyError(f"No material specified for section {name}.")

            if name in model.sections:
                raise KeyError("Redundant definition for section f{name}")

            Section = getSectionClass(sectionType)

            theSection = Section(name, data, materialName, model, **definition)

            model.sections[name] = theSection

        return model

    def createAnalyticalFieldsFromInputFile(self, model, inputFile):
        """Collects field defintions from the input file.

        Parameters
        ----------
        model
            A dictionary containing the model tree.
        inputFile
            A dictionary contaning the input file tree.

        Returns
        -------
        dict
            The updated model tree.
        """

        for fieldDef in inputFile["*analyticalField"]:
            analyticalFieldName = fieldDef["name"]
            analyticalFieldType = fieldDef["type"]
            analyticalFieldData = fieldDef["data"]

            analyticalFieldClass = getAnalyticalFieldByName(analyticalFieldType)
            analyticalField = analyticalFieldClass(analyticalFieldName, analyticalFieldData, model)

            model.analyticalFields[analyticalFieldName] = analyticalField

        return model
