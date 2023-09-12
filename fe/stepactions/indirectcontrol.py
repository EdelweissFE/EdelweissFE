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

"""
Indirect (displacement) controller for the NISTArcLength solver
"""

documentation = {
    "dof1": "Degree of freedom for the constraint ( model access expression )",
    "dof2": "Degree of freedom for the constraint ( model access expression )",
    "L": "Final distance (e.g. crack opening)",
}

from fe.stepactions.base.stepactionbase import StepActionBase
import numpy as np
from fe.utils.math import evalModelAccessibleExpression


class StepAction(StepActionBase):
    identification = "IndirectControl"

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name
        self.journal = journal
        self.model = model
        self.currentL0 = 0.0

        self.updateStepAction(action, jobInfo, model, fieldOutputController, journal)

    def computeDDLambda(self, dU, ddU_0, ddU_f, increment, dofManager):
        idcs = np.hstack([dofManager.idcsInDofVector[self.dof1], dofManager.idcsInDofVector[self.dof2]])

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        dL = incrementSize * self.L

        denom = self.c.dot(ddU_f[idcs])

        ddLambda = (dL - self.c.dot(dU[idcs] + ddU_0[idcs])) / self.c.dot(ddU_f[idcs])
        return ddLambda

    def finishIncrement(self, U, dU, dLambda, increment, dofManager):
        self.journal.message(
            "Dof 1: {:}, Dof 2: {:}".format(self.dof1.values, self.dof2.values),
            self.identification,
        )

    def applyAtStepEnd(self, model):
        self.currentL0 = self.c1.dot(self.dof1.values) + self.c2.dot(self.dof2.values)

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
