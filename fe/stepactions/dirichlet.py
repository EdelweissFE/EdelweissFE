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
    "nSet": "The node set for application of the BC",
    "1,2,3,...": "Prescribed values for components of the physical field",
    "components": "Prescribed values using a np.ndarray for representation; use 'x' for ignored values",
    "field": "Field for which the boundary condition is active",
    "analyticalField": "(Optional) scales the defined boundary condition",
    "f(t)": "(Optional) define an amplitude in the step progress interval [0...1]",
}


from fe.stepactions.base.stepactionbase import StepActionBase
from fe.config.phenomena import getFieldSize
import numpy as np
import sympy as sp


class StepAction(StepActionBase):
    """Dirichlet boundary condition, based on a node set"""

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):

        self.name = name

        dirichletIndices = []
        dirichletDelta = []

        nodeSets = model["nodeSets"]
        self.field = action["field"]

        if "analyticalField" in action:
            self.analyticalField = model["analyticalFields"][action["analyticalField"]]

        self.action = action
        self.nSet = nodeSets[action["nSet"]]

        self.possibleComponents = [str(i + 1) for i in range(getFieldSize(self.field, model["domainSize"]))]

        if "components" in action:
            action = self._getDirectionsFromComponents(action)

        for x, direction in enumerate(self.possibleComponents):
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

        if not dirichletIndices:
            raise ValueError("Invalid dirichlet components specified")

        self.indices = np.array(dirichletIndices)
        self.delta = np.array(dirichletDelta)

        self.amplitude = self._getAmplitude(action)

        self.active = True

    def applyAtStepEnd(self, U, P):

        self.active = False

    def updateStepAction(self, name, action, jobInfo, model, fieldOutputController, journal):

        self.active = True
        dirichletIndices = []
        dirichletDelta = []

        if "components" in action:
            action = self._getDirectionsFromComponents(action)

        for x, direction in enumerate(self.possibleComponents):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                dirichletIndices += directionIndices
                dirichletDelta += [float(action[direction])] * len(directionIndices)

        self.indices = np.array(dirichletIndices, dtype=int)

        self.delta = np.array(dirichletDelta)

        self.amplitude = self._getAmplitude(action)

    def getDelta(self, increment):

        if self.active:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            return self.delta * (self.amplitude(stepProgress) - (self.amplitude(stepProgress - incrementSize)))
        else:
            return 0.0

    def _getDirectionsFromComponents(self, action: dict) -> dict:
        """Determine the direction components from a numpy array representation.

        Parameters
        ----------
        action
            The dictionary defining this step action.

        Returns
        -------
        dict
            The updated dictionary defining this step action containing the directional definitions.
        """

        components = np.array(eval(action["components"].replace("x", "np.nan")), dtype=float)

        for i, t in enumerate(components):
            if not np.isnan(t):
                action[str(i + 1)] = t

        return action

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
