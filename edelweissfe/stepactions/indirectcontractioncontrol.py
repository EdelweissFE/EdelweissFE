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
# Created on Thu May 12 18:35:44 2022

# @author: Matthias Neuner

import numpy as np

from edelweissfe.stepactions.base.stepactionbase import StepActionBase
from edelweissfe.timesteppers.timestep import TimeStep

"""
Indirect (displacement) controller for the NISTArcLength solver
uses a ring to control the contraction, e.g., for tunneling simulations.

Currently 2D only!

The center is autotically computed from the bounding node coordinates.
"""

documentation = {
    "contractionNSet": "The node set defining the contraction ring",
    "L": "Final distance (e.g. crack opening)",
    "exportCVector": "(Optional) file to export the computed c vector",
}


class StepAction(StepActionBase):
    identification = "IndirectControl"

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name
        self.journal = journal

        self.currentL0 = 0.0

        self.L = float(action["L"])

        self.generateCVectorAndIndices(name, action, jobInfo, model, fieldOutputController, journal)

        if "exportCVector" in action:
            np.savetxt(action["exportCVector"] + ".csv", self.cVector)

        self.definition = str(action.get("definition", "absolute"))

    def computeDDLambda(self, dU, ddU_0, ddU_f, timeStep: TimeStep):
        dL = timeStep.stepProgressIncrement * self.L

        ddLambda = (dL - self.cVector.dot(dU[self.idcs] + ddU_0[self.idcs])) / self.cVector.dot(ddU_f[self.idcs])
        return ddLambda

    def finishIncrement(self, U, dU, dLambda):
        pass

    def applyAtStepEnd(self, U, P):
        self.currentL0 = self.cVector.dot(U[self.idcs])

    def updateStepAction(self, action, jobInfo, model, fieldOutputController, journal):
        if self.definition == "absolute":
            self.L = float(action["L"]) - self.currentL0
        else:
            self.L = float(action["L"])

        self.generateCVectorAndIndices(action, jobInfo, model, fieldOutputController, journal)

    def generateCVectorAndIndices(self, action, jobInfo, model, fieldOutputController, journal):
        contractionNSet = model.nodeSets[action["contractionNSet"]]

        nNodes = len(contractionNSet)

        allCoordinates = np.array([n.coordinates for n in contractionNSet])

        x_min = np.min(allCoordinates[:, 0])
        x_max = np.max(allCoordinates[:, 0])
        y_min = np.min(allCoordinates[:, 1])
        y_max = np.max(allCoordinates[:, 1])

        x_center = 0.5 * (x_max + x_min)
        y_center = 0.5 * (y_max + y_min)

        cVector = []
        idcsInDofVector = []

        for n in contractionNSet:
            vec_n_to_center = np.array([x_center - n.coordinates[0], y_center - n.coordinates[1]])
            norm_vec_n_to_center = np.linalg.norm(vec_n_to_center)

            vec_n_to_center_normalized = vec_n_to_center / norm_vec_n_to_center

            cVector.append(vec_n_to_center_normalized)

            idcsInDofVector.append([n.fields["displacement"][dim] for dim in range(2)])

        self.cVector = np.hstack(cVector)

        # dividing c vector to make 'average' contraction of ring:
        self.cVector *= 1.0 / nNodes

        self.idcs = np.hstack(idcsInDofVector)
