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

from edelweissfe.stepactions.base.stepactionbase import StepActionBase
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.config.phenomena import getFieldSize
import numpy as np
import sympy as sp
from abc import ABC, abstractmethod
from edelweissfe.sets.nodeset import NodeSet


class NodalLoadBase(StepActionBase):
    @property
    @abstractmethod
    def field(self) -> str:
        pass

    @property
    @abstractmethod
    def nodeSet(self) -> NodeSet:
        pass

    @abstractmethod
    def getCurrentLoad(self, timeStep: TimeStep) -> np.ndarray:
        """Return the current magnitude for this body load.

        Returns
        -------
        np.ndarray
            The magnitude.
        """
        pass