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

from abc import ABC, abstractmethod

import numpy as np

from edelweissfe.models.femodel import FEModel


class AnalyticalField(ABC):
    @abstractmethod
    def __init__(self, name: str, data: list[str], model: FEModel):
        """The field base class.

        Parameters
        ----------
        name
            The name of the field.
        data
            The input file data lines.
        model
            The dictionary containing the model tree.
        """

    @abstractmethod
    def evaluateAtCoordinates(self, coordinates: np.ndarray):
        """Evaluate a field at coordinates.

        Parameters
        ----------
        coordinates
            The spatial coordinates for evaluation of the field.

        Returns
        -------
        float
            The value of the field.
        """
