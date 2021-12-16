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
Created on Tue Jan 24 19:33:06 2017

@author: Matthias Neuner

Apply simple node forces on a nSet.
"""

documentation = {
    "nSet": "nSet for application of the BC",
    "1,2,3": "prescribed values in directions",
    "field": "field for BC",
    "f(t)": "(optional) define an amplitude",
}

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np
import sympy as sp


class StepAction(StepActionBase):
    """Defines node based load, defined on a nodeset."""

    def __init__(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.name = name
        nodeForceIndices = []
        nodeForceDelta = []
        nodeSets = modelInfo["nodeSets"]

        self.field = action["field"]
        self.idle = False
        self.nSet = nodeSets[action["nSet"]]

        for x, direction in enumerate(["1", "2", "3"]):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                nodeForceIndices += directionIndices
                nodeForceDelta += [float(action[direction])] * len(directionIndices)

        self.indices = np.asarray(nodeForceIndices, dtype=np.int)
        self.nodeForcesStepStart = np.zeros_like(self.indices, dtype=np.double)
        self.nodeForcesDelta = np.asarray(nodeForceDelta)
        self.currentNodeForces = np.zeros_like(self.nodeForcesDelta)

        self.amplitude = self.getAmplitude(action)

    def getAmplitude(self, action):

        if "f(t)" in action:
            t = sp.symbols("t")
            amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            amplitude = lambda x: x

        return amplitude

    def finishStep(self, U, P, stepMagnitude=None):

        if not self.idle:
            if stepMagnitude == None:
                # standard case
                self.nodeForcesStepStart += self.nodeForcesDelta * self.amplitude(1.0)
            else:
                # set the 'actual' increment manually, e.g. for arc length method
                self.nodeForcesStepStart += self.nodeForcesDelta * stepMagnitude

            self.nodeForcesDelta = 0
            self.idle = True

    def updateStepAction(self, action):

        self.idle = False
        nodeForceDelta = []
        for x, direction in enumerate(["1", "2", "3"]):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                nodeForceDelta += [float(action[direction])] * len(directionIndices)

        self.nodeForcesDelta = np.asarray(nodeForceDelta)
        self.amplitude = self.getAmplitude(action)

    def applyOnP(self, P, increment):

        if self.idle:
            P[self.indices] += self.nodeForcesStepStart
        else:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            t = stepProgress
            amp = self.amplitude(t)
            P[self.indices] += self.nodeForcesStepStart + self.nodeForcesDelta * amp

        return P
