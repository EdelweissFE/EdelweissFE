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
Stepaction to change material properties.

"""

documentation = {
    "material": "the id of the material to be changed",
    "index": "the index of the property in the material properties vector",
    "f(t)": "define a function depending on time"
}

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np
import sympy as sp

class StepAction(StepActionBase):
    """Action class for changing a material property"""

    def __init__(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.name = name
        self.theIndex = int (action ['index'])
        self.theMaterial = modelInfo['materials'][action['material']]

        self.updateStepAction(name, action, jobInfo, modelInfo, fieldOutputController, journal)

    def finishStep(self, U, P, stepMagnitude=None):
        self.active = False
        pass 

    def updateStepAction(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):
        """Update the function describing the material property"""
        self.active = True

        if "f(t)" in action:
            t = sp.symbols("t")
            self.f_t= sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")

    def changeTheProperties(self, increment, journal):
        """Change the actual properties depending on the current step time"""

        incNumber, incrementSize, stepProgress, dT, stepTime_n, totalTime_n = increment

        theCurrentProperty = self.f_t(stepTime_n + dT)
        journal.message(
                "Changing property[{:}] of material {:} to {:}".format(self.theIndex, self.theMaterial["name"], theCurrentProperty),
            self.name,
        )

        self.theMaterial['properties'][self.theIndex] = theCurrentProperty

