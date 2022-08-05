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
# Created on Sun Sep 10 14:10:20 2017

# @author: Matthias Neuner
"""
Module for applying displacements on a shotcrete shell
"""

documentation = {
    "nSet": "nSet for application of the BC",
    "displacements": "file containing the node displacements over time (column 0)",
}

from fe.stepactions.base.stepactionbase import StepActionBase
import numpy as np


class StepAction(StepActionBase):
    """Dirichlet boundary condition, based on a node set"""

    def __init__(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.name = name

        dirichletIndices = []

        nodeSets = modelInfo["nodeSets"]

        displacementsFile = action["displacements"]
        nSet = action["nSet"]

        x = np.loadtxt(displacementsFile)

        self.t = x[:, 0]
        self.t[:] -= self.t[0]

        self.U = x[:, 1:]

        nodes = nodeSets[nSet]

        dirichletIndices = [node.fields["displacement"] for node in nodes]
        self.indices = np.array(dirichletIndices).ravel()

    def finishStep(
        self,
    ):
        pass

    def updateStepAction(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):
        pass

    def getDelta(self, increment):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        if dT == 0.0:
            return np.zeros_like(self.indices)

        delta = np.array([np.interp((totalTime + dT), self.t, x) - np.interp((totalTime), self.t, x) for x in self.U.T])

        return delta
