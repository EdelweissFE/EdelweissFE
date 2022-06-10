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
#  Paul Hofer Paul.Hofer@uibk.ac.at
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

documentation = {
    "f(x,y,z)": "python expression using variables x, y, z (coordinates); dictionaries contained in modelInfo can be accessed",
}

from fe.utils.misc import stringDict
from fe.analyticalFields.base.analyticalFieldBase import (
    AnalyticalField as AnalyticalFieldBase,
)
from fe.utils.math import createModelAccessibleFunction


class AnalyticalField(AnalyticalFieldBase):
    def __init__(self, name, data, modelInfo):
        self.name = name
        self.type = "scalarExpression"

        self.domainSize = modelInfo["domainSize"]
        self.options = stringDict([e for l in data for e in l])

        expressionString = self.options["f(x,y,z)"]

        self.expression = createModelAccessibleFunction(
            expressionString, modelInfo, *"xyz"[: self.domainSize]
        )

        return

    def evaluateAtCoordinates(self, coords):
        for i1 in range(self.domainSize - len(coords)):
            coords.append(0)
        return self.expression(*coords)
