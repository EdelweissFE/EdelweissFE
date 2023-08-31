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
"""Define a field using a sclar expression.
"""

documentation = {
    "f(x,y,z)": "Python expression using variables x, y, z (coordinates); dictionaries contained in model can be accessed",
}

from fe.utils.misc import convertLinesToStringDictionary
from fe.analyticalfields.base.analyticalfieldbase import (
    AnalyticalField as AnalyticalFieldBase,
)
from fe.utils.math import createModelAccessibleFunction


class AnalyticalField(AnalyticalFieldBase):
    def __init__(self, name, data, model):
        self.name = name
        self.type = "scalarExpression"

        self.domainSize = model.domainSize
        self.options = convertLinesToStringDictionary(data)

        expressionString = self.options["f(x,y,z)"]

        self.expression = createModelAccessibleFunction(expressionString, model, *"xyz"[: self.domainSize])

        return

    def evaluateAtCoordinates(self, coords):
        for i1 in range(self.domainSize - len(coords)):
            coords.append(0)
        return self.expression(*coords)
