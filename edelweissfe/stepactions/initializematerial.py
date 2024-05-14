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
# Created on Mo July 29 10:50:53 2019

# @author: Matthias Neuner
"""
Let materials initialize themselves (e.g., state vars depending on material parameters...) !
"""

import numpy as np

from edelweissfe.stepactions.base.stepactionbase import StepActionBase

documentation = {}


class StepAction(StepActionBase):
    """Initializes materials"""

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name

        self.theElements = model.elementSets[action.get("elSet", "all")]
        self.active = True
        self.emptyDef = np.array([0.0])

    def applyAtStepEnd(self, model, stepMagnitude=None):
        self.active = False

    def applyAtStepStart(self, model):
        if not self.active:
            return

        for el in self.theElements:
            el.setInitialCondition("initialize material", self.emptyDef)
