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
#  Paul Hofer paul.hofer@uibk.ac.at
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
#
# @author: Matthias Neuner, Paul Hofer

from abc import ABC, abstractmethod

import numpy as np

from edelweissfe.utils.math import createFunction
from edelweissfe.utils.misc import (
    convertAssignmentsToStringDictionary,
    splitLineAtCommas,
    strCaseCmp,
)


class Section(ABC):
    def __init__(self, name, datalines, materialName, model, **kwargs):
        self.materialParameterFromFieldDefs = []
        self.writeMaterialPropertiesToFile = False
        elSetNames = []

        for line in datalines:
            line = splitLineAtCommas(line)

            if strCaseCmp(line[0], "materialParameterFromField"):
                definition = convertAssignmentsToStringDictionary(line[1:])

                if "type" in definition.keys():
                    if not any(
                        [strCaseCmp(definition["type"], implementedType) for implementedType in ["setToValue", "scale"]]
                    ):
                        raise KeyError(
                            "{}: {} is not a known type; currently available types: 'setToValue', 'scale'".format(
                                name, definition["type"]
                            )
                        )
                else:
                    raise KeyError(
                        "{}: Type option must be set; currently available types: 'setToValue', 'scale'".format(name)
                    )

                if "f(p,f)" in definition.keys():
                    definition["expression"] = createFunction(definition["f(p,f)"], "p", "f", model=model)
                else:
                    definition["expression"] = createFunction("f", "p", "f", model=model)

                self.materialParameterFromFieldDefs.append(definition)
            elif strCaseCmp(line[0], "writeMaterialPropertiesToFile"):
                self.writeMaterialPropertiesToFile = True
                definition = convertAssignmentsToStringDictionary(line[1:])
                self.materialPropertiesFileName = definition.get("filename")
            else:
                elSetNames.extend(line)

        self.elSets = [model.elementSets[setName] for setName in elSetNames]
        self.material = model.materials[materialName]

    def assignSectionPropertiesToModel(self, model):
        if any(self.materialParameterFromFieldDefs):
            for elSet in self.elSets:
                for el in elSet:
                    modifiedMaterial = self.material.copy()
                    modifiedMaterial["properties"] = self.propertiesFromField(el, self.material, model)

                    self.assignSectionPropertiesToElement(el, material=modifiedMaterial)

        else:
            for elSet in self.elSets:
                for el in elSet:
                    self.assignSectionPropertiesToElement(el)

        if self.writeMaterialPropertiesToFile:
            self.exportMaterialPropertiesToFile(self.elSets)

        return model

    @abstractmethod
    def assignSectionPropertiesToElement(self, element, **kwargs):
        pass

    def propertiesFromField(self, el, material, model):
        coordinatesAtCenter = el.getCoordinatesAtCenter()
        materialProperties = np.copy(material["properties"])

        for definition in self.materialParameterFromFieldDefs:
            index = int(definition["index"])
            fieldValue = model.analyticalFields[definition["field"]].evaluateAtCoordinates(coordinatesAtCenter)[0][0]
            parameterValue = materialProperties[index]

            if strCaseCmp(definition["type"], "setToValue"):
                materialProperties[index] = definition["expression"](parameterValue, fieldValue)
            elif strCaseCmp(definition["type"], "scale"):
                materialProperties[index] *= definition["expression"](parameterValue, fieldValue)

        return materialProperties

    def exportMaterialPropertiesToFile(self, elSets):
        with open("{:}.csv".format(self.materialPropertiesFileName), "w+") as f:
            for elSet in elSets:
                for el in elSet:
                    f.write("{:}".format(el.elNumber))
                    [f.write("{:} ".format(matprop)) for matprop in el._materialProperties]
                    f.write("\n")
