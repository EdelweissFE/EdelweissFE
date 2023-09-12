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
# Created on Tue Jan 24 19:33:06 2017

# @author: Matthias Neuner

"""
Apply simple node forces on a nSet.
"""

documentation = {
    "nSet": "Node set for application of the boundary condition",
    "1,2,3,...": "Prescribed values for components of the physical field",
    "components": "Prescribed values using a np.ndarray for representation; use 'x' for ignored values",
    "field": "Field for which the boundary condition is active",
    "f(t)": "(Optional) define an amplitude in the step progress interval [0...1]",
}

from fe.stepactions.base.stepactionbase import StepActionBase
from fe.config.phenomena import getFieldSize
import numpy as np
import sympy as sp


class StepAction(StepActionBase):
    """Defines node based load, defined on a nodeset."""

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name
        # nodeForceIndices = []
        # nodeForceDelta = []
        nodeSets = model.nodeSets

        self.field = action["field"]
        # self.idle = False
        self.nSet = nodeSets[action["nSet"]]

        self.fieldSize = getFieldSize(self.field, model.domainSize)

        shape = (len(self.nSet), self.fieldSize)

        self.nodeForcesStepStart = np.zeros(shape)
        self.nodeForcesDelta = np.zeros(shape)

        self.possibleComponents = [str(i + 1) for i in range(self.fieldSize)]

        self.updateStepAction(action, jobInfo, model, fieldOutputController, journal)

        # if "components" in action:
        #     action = self._getDirectionsFromComponents(action)

        # for x, direction in enumerate(self.possibleComponents):
        #     directionIndices = [node.fields[self.field][x] for node in self.nSet]
        #     if direction in action:
        #         nodeForceIndices += directionIndices
        #         nodeForceDelta += [float(action[direction])] * len(directionIndices)

        # self.indices = np.asarray(nodeForceIndices, dtype=int)
        # self.nodeForcesStepStart = np.zeros_like(self.indices, dtype=np.double)
        # self.nodeForcesDelta = np.asarray(nodeForceDelta)
        # self.currentNodeForces = np.zeros_like(self.nodeForcesDelta)

    def updateStepAction(self, action, jobInfo, model, fieldOutputController, journal):
        """Update the step action.

        It is a reasonable requirement that the updated direction components cannot change.
        """

        self.idle = False
        # nodeForceDelta = []
        # nodeForceIndices = []

        if "components" in action:
            nodeLoad = np.asarray(eval(action["components"].replace("x", "0")), dtype=float)
        else:
            nodeLoad = self._getComponentsFromDirection(action)

        # for x, direction in enumerate(self.possibleComponents):
        #     # directionIndices = [node.fields[self.field][x] for node in self.nSet]
        #     if direction in action:
        #         # nodeForceIndices += directionIndices
        #         nodeForceDelta += [float(action[direction])] * len(directionIndices)

        # if not (nodeForceIndices == self.indices).all():
        #     raise ValueError("Components for node forces action can not change!")

        nodeForcesDelta = np.tile(nodeLoad, (len(self.nSet), 1))

        self.nodeForcesDelta = nodeForcesDelta
        self.amplitude = self._getAmplitude(action)

        # self.amplitude = self._getAmplitude(action)

    def _getAmplitude(self, action: dict) -> callable:
        """Determine the amplitude for the step, depending on a potentially specified function.

        Parameters
        ----------
        action
            The dictionary defining this step action.

        Returns
        -------
        callable
            The function defining the amplitude depending on the step propress.
        """

        if "f(t)" in action:
            t = sp.symbols("t")
            amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            amplitude = lambda x: x

        return amplitude

    def applyAtStepEnd(self, model, stepMagnitude=None):
        if not self.idle:
            if stepMagnitude == None:
                # standard case
                self.nodeForcesStepStart += self.nodeForcesDelta * self.amplitude(1.0)
            else:
                # set the 'actual' increment manually, e.g. for arc length method
                self.nodeForcesStepStart += self.nodeForcesDelta * stepMagnitude

            self.nodeForcesDelta[:] = 0
            self.idle = True

    def getLoad(self, increment):
        if self.idle:
            return self.nodeForcesStepStart
            # P[self.indices] += self.nodeForcesStepStart
        else:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            t = stepProgress
            amp = self.amplitude(t)

            return self.nodeForcesStepStart + self.nodeForcesDelta * amp

    def _getComponentsFromDirection(self, action: dict) -> np.ndarray:
        nodeLoad = np.zeros(self.fieldSize)

        for i, comp in enumerate(self.possibleComponents):
            if comp in action:
                nodeLoad[i] = float(action[comp])

        return nodeLoad
