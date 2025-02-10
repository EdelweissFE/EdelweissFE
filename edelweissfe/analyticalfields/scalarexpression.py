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
"""Define a field using a sclar expression."""

import numpy as np

from edelweissfe.analyticalfields.base.analyticalfieldbase import (
    AnalyticalField as AnalyticalFieldBase,
)
from edelweissfe.utils.math import createModelAccessibleFunction
from edelweissfe.utils.misc import convertLinesToStringDictionary

documentation = {
    "f(x,y,z)": "Python expression using variables x, y, z (coordinates); dictionaries contained in model can be accessed",
}


class AnalyticalField(AnalyticalFieldBase):
    def __init__(self, name, data, model):
        self.name = name
        self.type = "scalarExpression"

        self.domainSize = model.domainSize
        self.options = convertLinesToStringDictionary(data)

        expressionString = self.options["f(x,y,z)"]

        self.expression = createModelAccessibleFunction(expressionString, model, *"xyz")  # [: self.domainSize])

        return

    def evaluateAtCoordinates(self, coords):
        coords = np.array(coords)

        if coords.ndim == 1:
            coords = np.expand_dims(coords, 0)
        coords = np.c_[coords, np.zeros((coords.shape[0], 3 - coords.shape[-1]))]

        return np.expand_dims(np.array([float(self.expression(*coords_)) for coords_ in coords]), 1)
