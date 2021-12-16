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
"""
Created on Tue Feb 9 10:05:41 2021

@author: Matthias Neuner
"""

documentation = {
    "fieldOutput": "field output to set",
    "type": "const",
    "value": "values",
}

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np


class StepAction(StepActionBase):
    """Defines node based load, defined on a nodeset."""

    def __init__(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.name = name
        self.active = True
        self.journal = journal
        self.fieldOutputName = action["fieldOutput"]
        self.fieldOutput = fieldOutputController.fieldOutputs[self.fieldOutputName]
        self.type = action["type"]
        self.value = action["value"]

    def finishStep(self, U, P, stepMagnitude=None):
        self.active = False

    def updateStepAction(self, action):
        self.active = True

    def apply(
        self,
    ):
        if self.type == "uniform":
            currentResults = np.zeros_like(self.fieldOutput.getLastResult())
            currentResults[:] = float(self.value)
            self.journal.message(
                "setting field {:} to uniform value {:}".format(self.fieldOutputName, self.value), self.identification
            )
            self.fieldOutput.setResults(currentResults)
