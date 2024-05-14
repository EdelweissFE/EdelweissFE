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

import numpy as np
import sympy as sp

from edelweissfe.config.phenomena import getFieldSize
from edelweissfe.stepactions.base.dirichletbase import DirichletBase
from edelweissfe.timesteppers.timestep import TimeStep

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


class StepAction(DirichletBase):
    """Dirichlet boundary condition, based on a node set"""

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name

        self.field = action["field"]

        self.action = action
        self.nSet = model.nodeSets[action["nSet"]]
        self.fieldSize = getFieldSize(self.field, model.domainSize)
        self.possibleComponents = [str(i + 1) for i in range(self.fieldSize)]

        self._components = None

        self.updateStepAction(action, jobInfo, model, fieldOutputController, journal)

    @property
    def components(
        self,
    ):
        return self._components

    def applyAtStepEnd(self, model):
        self.active = False

    def updateStepAction(self, action, jobInfo, model, fieldOutputController, journal):
        self.active = True

        self.action = action

        if "components" in action:
            action = self._getDirectionsFromComponents(action)

        components = [i for i, direction in enumerate(self.possibleComponents) if direction in action]

        values = {
            i: float(action[direction]) for i, direction in enumerate(self.possibleComponents) if direction in action
        }

        self.delta = np.tile(list(values.values()), (len(self.nSet), 1))

        # for i, node in enumerate(self.nSet):
        if "analyticalField" in action:
            self.analyticalField = model.analyticalFields[action["analyticalField"]]
            for i, node in enumerate(self.nSet):
                self.delta[i, :] *= self.analyticalField.evaluateAtCoordinates(node.coordinates)[0][0]

        self._components = components

        self.amplitude = self._getAmplitude(action)

    def getDelta(self, timeStep: TimeStep):
        if self.active:
            return self.delta * (
                self.amplitude(timeStep.stepProgress)
                - (self.amplitude(timeStep.stepProgress - timeStep.stepProgressIncrement))
            )
        else:
            return self.delta * 0.0

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

            def amplitude(x):
                return x

        return amplitude
