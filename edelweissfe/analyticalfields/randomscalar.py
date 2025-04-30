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
"""Define a random field using the GSTools library."""

import gstools
import numpy as np

from edelweissfe.analyticalfields.base.analyticalfieldbase import (
    AnalyticalField as AnalyticalFieldBase,
)
from edelweissfe.utils.caseinsensitivedict import CaseInsensitiveDict
from edelweissfe.utils.misc import convertLinesToStringDictionary, strCaseCmp

documentation = {
    "model": "(Optional) default = Gaussian",
    "mean": "(Optional) default = 0.",
    "variance": "(Optional) default = 1.",
    "lengthScale": "(Optional) default = 10.",
    "seed": "(Optional) default = 0",
}


class AnalyticalField(AnalyticalFieldBase):
    def __init__(self, name, data, model):
        self.name = name
        self.type = "randomScalar"

        options = CaseInsensitiveDict(convertLinesToStringDictionary(data))

        self.domainSize = model.domainSize

        modelType = options.get("model", "Gaussian")

        if strCaseCmp(modelType, "Gaussian"):
            variance = float(options.get("variance", 1.0))
            lengthScale = float(options.get("lengthScale", 10.0))

            # modelMethod = getattr(gstools, modelType)
            model = gstools.Gaussian(
                dim=self.domainSize,
                var=variance,
                len_scale=lengthScale,
            )
        elif strCaseCmp(modelType, "Matern"):
            mean = float(options.get("mean", 0.0))
            variance = float(options.get("variance", 1.0))
            lengthScale = float(options.get("lengthScale", 10.0))
            nu = float(options.get("nu", 1.0))

            # modelMethod = getattr(gstools, modelType)
            model = gstools.covmodel.Matern(
                dim=self.domainSize,
                var=variance,
                len_scale=lengthScale,
                nu=nu,
            )
        else:
            raise NotImplementedError(f"Model type {modelType} not implemented.")

        mean = float(options.get("mean", 0.0))
        seed = int(options.get("seed", 0))
        self.srf = gstools.SRF(model, seed=seed, mean=mean)

        return

    def evaluateAtCoordinates(self, coords):
        coords = np.array(coords)

        if coords.ndim == 1:
            coords = np.expand_dims(coords, 0)

        return np.expand_dims(np.array([self.srf(coords_)[0] for coords_ in coords]), 1)
