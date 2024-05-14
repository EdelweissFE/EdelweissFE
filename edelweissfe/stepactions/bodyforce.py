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

import numpy as np
import sympy as sp

from edelweissfe.stepactions.base.bodyloadbase import BodyLoadBase
from edelweissfe.timesteppers.timestep import TimeStep

"""
Simple body force load.
If not modified in subsequent steps, the load held constant.
"""

documentation = {
    "forceVector": "The force vector",
    "delta": "In subsequent steps only: define the updated force vector incrementally",
    "f(t)": "(Optional) define an amplitude in the step progress interval [0...1]",
}


class StepAction(BodyLoadBase):
    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self._name = name
        self._forceAtStepStart = 0.0
        self._elSet = model.elementSets[action["elSet"]]
        load = np.fromstring(action["forceVector"], sep=",", dtype=np.double)

        if len(load) < model.domainSize:
            raise Exception("BodyForce {:}: force vector has wrong dimension!".format(self._name))

        self._delta = load
        if "f(t)" in action:
            t = sp.symbols("t")
            self._amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            self._amplitude = lambda x: x

        self._idle = False

    def applyAtStepEnd(self, model, stepMagnitude=None):
        if not self._idle:
            if stepMagnitude is None:
                # standard case
                self._forceAtStepStart += self._delta * self._amplitude(1.0)
            else:
                # set the 'actual' increment manually, e.g. for arc length method
                self._forceAtStepStart += self._delta * stepMagnitude

            self._delta = 0
            self._idle = True

    def updateStepAction(self, action, jobInfo, model, fieldOutputController, journal):
        if "forceVector" in action:
            self._delta = np.fromstring(action["forceVector"], sep=",", dtype=np.double) - self._forceAtStepStart
        elif "delta" in action:
            self._delta = np.fromstring(action["delta"], sep=",", dtype=np.double)

        if "f(t)" in action:
            t = sp.symbols("t")
            self._amplitude = sp.lambdify(t, sp.sympify(action["f(t)"]), "numpy")
        else:
            self._amplitude = lambda x: x

        self._idle = False

    def getCurrentLoad(self, timeStep: TimeStep):
        if self._idle is True:
            t = 1.0
        else:
            t = timeStep.stepProgress

        return self._forceAtStepStart + self._delta * self._amplitude(t)

    @property
    def elementSet(self) -> list:
        return self._elSet
