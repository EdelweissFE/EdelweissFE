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
# Created on Thu Nov  2 18:35:44 2017

# @author: Matthias Neuner

import numpy as np

from edelweissfe.numerics.dofmanager import DofManager
from edelweissfe.stepactions.base.stepactionbase import StepActionBase
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.utils.math import evalModelAccessibleExpression

"""
Indirect (displacement) controller for the NISTArcLength solver
"""

documentation = {
    "dof1": "Degree of freedom for the constraint ( model access expression )",
    "dof2": "Degree of freedom for the constraint ( model access expression )",
    "L": "Final distance (e.g. crack opening)",
}


class StepAction(StepActionBase):
    identification = "IndirectControl"

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name
        self.journal = journal
        self.model = model
        self.currentL0 = 0.0

        self.updateStepAction(action, jobInfo, model, fieldOutputController, journal)

    def computeDDLambda(self, dU, ddU_0, ddU_f, timeStep: TimeStep, dofManager: DofManager):
        idcs = np.hstack(
            [
                dofManager.idcsOfFieldVariablesInDofVector[self.dof1],
                dofManager.idcsOfFieldVariablesInDofVector[self.dof2],
            ]
        )

        dL = timeStep.stepProgressIncrement * self.L

        ddLambda = (dL - self.c.dot(dU[idcs] + ddU_0[idcs])) / self.c.dot(ddU_f[idcs])
        return ddLambda

    def finishIncrement(self, U, dU, dLambda, timeStep: TimeStep, dofManager):
        self.journal.message(
            f"C1·DOF1: {self.c1.dot(self.dof1.values)}, C2·DOF2: {self.c2.dot(self.dof2.values)}",
            self.identification,
        )

    def applyAtStepEnd(self, model):
        # self.currentL0 = self.c1.dot(self.dof1.values) + self.c2.dot(self.dof2.values)
        self.currentL0 = self.L

    def updateStepAction(self, action, jobInfo, model, fieldOutputController, journal):
        self.definition = str(action.get("definition", "absolute"))

        if self.definition == "absolute":
            self.L = float(action["L"]) - self.currentL0
        else:
            self.L = float(action["L"])

        self.dof1 = evalModelAccessibleExpression(action["dof1"], model)
        self.dof2 = evalModelAccessibleExpression(action["dof2"], model)

        self.c1 = np.asarray(eval(action["cVector1"].replace("x", "0")), dtype=float)
        self.c2 = np.asarray(eval(action["cVector2"].replace("x", "0")), dtype=float)

        self.c = np.hstack([self.c1, self.c2])
