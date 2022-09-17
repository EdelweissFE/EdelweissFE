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

import numpy as np
from fe.utils.misc import convertAssignmentsToStringDictionary, splitLineAtCommas
from fe.utils.misc import strCaseCmp
from fe.utils.math import createFunction

from abc import ABC, abstractmethod


class Section(ABC):
    def __init__(self, name, datalines, materialName, thickness, model):
        self.elSetNames = []
        self.materialParameterFromFieldDefs = []

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
            else:
                self.elSetNames.extend(line)

        self.materialName = materialName

    @abstractmethod
    def assignSectionPropertiesToModel(self, model):
        pass

    def propertiesFromField(self, el, material, model):
        coordinatesAtCenter = el.getCoordinatesAtCenter()
        materialProperties = np.copy(material["properties"])

        for definition in self.materialParameterFromFieldDefs:
            index = int(definition["index"])
            fieldValue = model["analyticalFields"][definition["field"]].evaluateAtCoordinates(coordinatesAtCenter)
            parameterValue = materialProperties[index]

            if strCaseCmp(definition["type"], "setToValue"):
                materialProperties[index] = definition["expression"](parameterValue, fieldValue)
            elif strCaseCmp(definition["type"], "scale"):
                materialProperties[index] *= definition["expression"](parameterValue, fieldValue)

        return materialProperties
