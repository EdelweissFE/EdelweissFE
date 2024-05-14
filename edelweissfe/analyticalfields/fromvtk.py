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
"""Use PyVista to interpolate from vtk data.
"""

import numpy as np
import pyvista

from edelweissfe.analyticalfields.base.analyticalfieldbase import (
    AnalyticalField as AnalyticalFieldBase,
)
from edelweissfe.utils.misc import convertLinesToStringDictionary

documentation = {
    "file": "path to database file",
    "result": "result name in database (optional if database contains only one dataset)",
}

# from inspect import signature


class AnalyticalField(AnalyticalFieldBase):
    """ """

    def __init__(self, name, data, model):
        self.name = name
        self.type = "fromVtk"

        options = convertLinesToStringDictionary(data)

        self.domainSize = model.domainSize

        file = options["file"]

        reader = pyvista.get_reader(file)
        self.data = reader.read()

        result = options.get("result")
        availableResults = self.data.array_names

        try:
            assert len(availableResults) > 0
        except AssertionError:
            raise AssertionError("Database does not contain at least one result.")

        if not result:
            if len(availableResults) == 1:
                result = self.data.array_names[0]
            else:
                raise ValueError("Database contains multiple results. Specify result with option result=...")
        else:
            if result not in availableResults:
                raise KeyError(f"Specified result '{result}' not available. Available results: {availableResults}")

        self.result = result

        return

    def evaluateAtCoordinates(self, coords):
        coords = np.array(coords)

        interpolatedData = pyvista.PointSet(coords).sample(self.data)
        interpolatedResult = interpolatedData[self.result]

        return np.expand_dims(interpolatedResult, 1)
