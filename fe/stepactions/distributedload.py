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
# Created on Wed May 10 13:12:40 2017

# @author: Matthias Neuner
"""
Standard distributed load, applied on a surface set.
If not modified in subsequent steps, the load held constant.
"""

documentation = {
    "surface": "Surface for application of the distributed load",
    "magnitude": "Magnitude of the distributed load",
    "delta": "In subsequent steps only: define the new magnitude incrementally",
    "f(t)": "(Optional) define an amplitude in the step progress interval [0...1]",
    "type": "The load type, e.g., pressure or surface traction; Must be supported by the element type",
}

from fe.stepactions.base.stepactionbase import StepActionBase
import numpy as np
import sympy as sp


class StepAction(StepActionBase):
    """Distributed load, defined on an element-based surface"""

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name
        self.magnitudeAtStepStart = 0.0
        self.surface = model.surfaces[action["surface"]]
        self.loadType = action["type"]
        magnitude = np.fromstring(action["magnitude"], sep=",")

        self.delta = magnitude
        if "f(t)" in action:
            t = sp.symbols("t")
            self.amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            self.amplitude = lambda x: x

        self.idle = False

    def applyAtStepEnd(self, model, stepMagnitude=None):
        if not self.idle:
            if stepMagnitude == None:
                # standard case
                self.magnitudeAtStepStart += self.delta * self.amplitude(1.0)
            else:
                # set the 'actual' increment manually, e.g. for arc length method
                self.magnitudeAtStepStart += self.delta * stepMagnitude

            self.delta = 0
            self.idle = True

    def updateStepAction(self, action, jobInfo, model, fieldOutputController, journal):
        if "magnitude" in action:
            self.delta = np.fromstring(action["magnitude"], sep=",") - self.magnitudeAtStepStart
        elif "delta" in action:
            self.delta = np.fromstring(action["delta"], sep=",")

        if "f(t)" in action:
            t = sp.symbols("t")
            self.amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            self.amplitude = lambda x: x

        self.idle = False

    def getCurrentMagnitude(self, increment):
        if self.idle == True:
            t = 1.0
        else:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            t = stepProgress

        return self.magnitudeAtStepStart + self.delta * self.amplitude(t)
