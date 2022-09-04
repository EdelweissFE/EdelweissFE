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
# Created on Thu Nov 15 13:15:14 2018

# @author: Matthias Neuner
"""
Simple body force load.
If not modified in subsequent steps, the load held constant.
"""

documentation = {
    "forceVector": "The force vector",
    "delta": "In subsequent steps only: define the updated force vector incrementally",
    "f(t)": "(Optional) define an amplitude in the step progress interval [0...1]",
}

from fe.stepactions.base.stepactionbase import StepActionBase
import numpy as np
import sympy as sp


class StepAction(StepActionBase):
    def __init__(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.name = name
        self.forceAtStepStart = 0.0
        self.elements = modelInfo["elementSets"][action["elSet"]]
        magnitude = np.fromstring(action["forceVector"], sep=",", dtype=np.double)

        if len(magnitude) < modelInfo["domainSize"]:
            raise Exception("BodyForce {:}: force vector has wrong dimension!".format(self.name))

        self.delta = magnitude
        if "f(t)" in action:
            t = sp.symbols("t")
            self.amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            self.amplitude = lambda x: x

        self.idle = False

    def applyAtStepEnd(self, U, P, stepMagnitude=None):

        if not self.idle:
            if stepMagnitude == None:
                # standard case
                self.forceAtStepStart += self.delta * self.amplitude(1.0)
            else:
                # set the 'actual' increment manually, e.g. for arc length method
                self.forceAtStepStart += self.delta * stepMagnitude

            self.delta = 0
            self.idle = True

    def updateStepAction(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        if "forceVector" in action:
            self.delta = np.fromstring(action["forceVector"], sep=",", dtype=np.double) - self.forceAtStepStart
        elif "delta" in action:
            self.delta = np.fromstring(action["delta"], sep=",", dtype=np.double)

        if "f(t)" in action:
            t = sp.symbols("t")
            self.amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            self.amplitude = lambda x: x

        self.idle = False

    def getCurrentBodyForce(self, increment):
        if self.idle == True:
            t = 1.0
        else:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            t = stepProgress

        return self.forceAtStepStart + self.delta * self.amplitude(t)
