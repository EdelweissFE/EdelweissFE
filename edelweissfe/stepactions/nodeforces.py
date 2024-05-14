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

import numpy as np
import sympy as sp

from edelweissfe.config.phenomena import getFieldSize
from edelweissfe.sets.nodeset import NodeSet
from edelweissfe.stepactions.base.nodalloadbase import NodalLoadBase
from edelweissfe.timesteppers.timestep import TimeStep

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


class StepAction(NodalLoadBase):
    """Defines node based load, defined on a nodeset."""

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name
        nodeSets = model.nodeSets

        self._field = action["field"]
        self._nSet = nodeSets[action["nSet"]]

        self._fieldSize = getFieldSize(self._field, model.domainSize)

        shape = (len(self._nSet), self._fieldSize)

        self.nodeForcesStepStart = np.zeros(shape)
        self.nodeForcesDelta = np.zeros(shape)

        self.possibleComponents = [str(i + 1) for i in range(self._fieldSize)]

        self.updateStepAction(action, jobInfo, model, fieldOutputController, journal)

    def updateStepAction(self, action, jobInfo, model, fieldOutputController, journal):
        """Update the step action.

        It is a reasonable requirement that the updated direction components cannot change.
        """

        self._idle = False

        if "components" in action:
            nodeLoad = np.asarray(eval(action["components"].replace("x", "0")), dtype=float)
        else:
            nodeLoad = self._getComponentsFromDirection(action)

        nodeForcesDelta = np.tile(nodeLoad, (len(self._nSet), 1))

        self.nodeForcesDelta = nodeForcesDelta
        self.amplitude = self._getAmplitude(action)

    @property
    def field(self) -> str:
        return self._field

    @property
    def nodeSet(self) -> NodeSet:
        return self._nSet

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

    def applyAtStepEnd(self, model, stepMagnitude=None):
        if not self._idle:
            if stepMagnitude is None:
                # standard case
                self.nodeForcesStepStart += self.nodeForcesDelta * self.amplitude(1.0)
            else:
                # set the 'actual' increment manually, e.g. for arc length method
                self.nodeForcesStepStart += self.nodeForcesDelta * stepMagnitude

            self.nodeForcesDelta[:] = 0
            self._idle = True

    def getCurrentLoad(self, timeStep: TimeStep):
        if self._idle:
            return self.nodeForcesStepStart
        else:
            t = timeStep.stepProgress
            amp = self.amplitude(t)

            return self.nodeForcesStepStart + self.nodeForcesDelta * amp

    def _getComponentsFromDirection(self, action: dict) -> np.ndarray:
        nodeLoad = np.zeros(self._fieldSize)

        for i, comp in enumerate(self.possibleComponents):
            if comp in action:
                nodeLoad[i] = float(action[comp])

        return nodeLoad
