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
"""Define a random field using the GSTools library.
"""

documentation = {
    "model": "(Optional) default = Gaussian",
    "mean": "(Optional) default = 0.",
    "variance": "(Optional) default = 1.",
    "lengthScale": "(Optional) default = 10.",
    "seed": "(Optional) default = 0",
}

import numpy as np
import sympy as sp
import gstools
from fe.utils.misc import convertLinesToStringDictionary
from fe.utils.math import createFunction
from fe.analyticalfields.base.analyticalfieldbase import (
    AnalyticalField as AnalyticalFieldBase,
)
from inspect import signature


class AnalyticalField(AnalyticalFieldBase):
    def __init__(self, name, data, model):
        self.name = name
        self.type = "randomScalar"

        options = convertLinesToStringDictionary(data)

        self.domainSize = model.domainSize

        modelType = options.get("model", "Gaussian")
        mean = float(options.get("mean", 0.0))
        variance = float(options.get("variance", 1.0))
        lengthScale = float(options.get("lengthScale", 10.0))
        seed = int(options.get("seed", 0))

        modelMethod = getattr(gstools, modelType)
        model = modelMethod(
            dim=self.domainSize,
            var=variance,
            len_scale=lengthScale,
        )
        self.srf = gstools.SRF(model, seed=seed, mean=mean)

        return

    def evaluateAtCoordinates(self, coords):
        for i1 in range(self.domainSize - len(coords)):
            coords.append(0)
        return float(self.srf(coords))
