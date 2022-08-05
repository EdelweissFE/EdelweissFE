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
# Created on Mon Jan 23 13:03:09 2017

# @author: Matthias Neuner
"""
Standard Dirichlet boundary condition.
If not modified in subsequent steps, the BC is held constant.
"""
documentation = {
    "nSet": "nSet for application of the BC",
    "1,2,3": "prescribed values in directions",
    "field": "field for BC",
    "analyticalField": "(optional) scales the defined BCs",
    "f(t)": "(optional) define an amplitude",
}

from fe.stepactions.base.stepactionbase import StepActionBase
import numpy as np
import sympy as sp


class StepAction(StepActionBase):
    """Dirichlet boundary condition, based on a node set"""

    def __init__(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.name = name

        dirichletIndices = []
        dirichletDelta = []

        nodeSets = modelInfo["nodeSets"]
        self.field = action["field"]

        if "analyticalField" in action:
            self.analyticalField = modelInfo["analyticalFields"][action["analyticalField"]]

        self.action = action
        self.nSet = nodeSets[action["nSet"]]

        for x, direction in enumerate(["1", "2", "3"]):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                dirichletIndices += directionIndices

                if "analyticalField" in action:
                    scaleFactors = [self.analyticalField.evaluateAtCoordinates(node.coordinates) for node in self.nSet]

                    unScaledDelta = [float(action[direction])] * len(directionIndices)
                    scaledDelta = []
                    for val1, val2 in zip(scaleFactors, unScaledDelta):
                        scaledDelta.append(val1 * val2)

                    dirichletDelta += scaledDelta

                else:
                    dirichletDelta += [float(action[direction])] * len(directionIndices)

        self.indices = np.array(dirichletIndices)
        self.delta = np.array(dirichletDelta)

        if "f(t)" in action:
            t = sp.symbols("t")
            self.amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            self.amplitude = lambda x: x

        self.active = True

    def finishStep(self, U, P):

        self.active = False

    def updateStepAction(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.active = True
        dirichletIndices = []
        dirichletDelta = []

        for x, direction in enumerate(["1", "2", "3"]):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                dirichletIndices += directionIndices
                dirichletDelta += [float(action[direction])] * len(directionIndices)

        self.delta = np.array(dirichletDelta)

        if "f(t)" in action:
            t = sp.symbols("t")
            self.amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            self.amplitude = lambda x: x

    def getDelta(self, increment):

        if self.active:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            return self.delta * (self.amplitude(stepProgress) - (self.amplitude(stepProgress - incrementSize)))
        else:
            return 0.0
