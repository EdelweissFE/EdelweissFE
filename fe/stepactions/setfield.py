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
# Created on Tue Feb 9 10:05:41 2021

# @author: Matthias Neuner
"""
Set a field (via fieldOutput) to a predefined value.
"""

documentation = {
    "fieldOutput": "field output to set",
    "type": "'const or 'analyticalField'",
    "value": "scalar value if type 'const'; name of analyticalField if type 'analyticalField'",
}

from fe.stepactions.base.stepactionbase import StepActionBase
import numpy as np


class StepAction(StepActionBase):
    def __init__(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.name = name
        self.active = True
        self.journal = journal
        self.fieldOutputName = action["fieldOutput"]
        self.fieldOutput = fieldOutputController.fieldOutputs[self.fieldOutputName]
        self.type = action["type"]
        self.value = action["value"]

        if not self.type in ["uniform", "analyticalField"]:
            raise Exception("Invalid type: {}".format(self.type))

        if self.type == "analyticalField":
            self.analyticalField = modelInfo["analyticalFields"][self.value]

    def finishStep(self, U, P, stepMagnitude=None):
        self.active = False

    def updateStepAction(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):
        self.active = True

    def apply(
        self,
    ):
        if self.type == "uniform":
            currentResults = np.zeros_like(self.fieldOutput.getLastResult())
            currentResults[:] = float(self.value)
            self.journal.message(
                "Setting field {:} to uniform value {:}".format(self.fieldOutputName, self.value), self.name
            )
            self.fieldOutput.setResults(currentResults)

        if self.type == "analyticalField":

            currentResults = np.zeros_like(self.fieldOutput.getLastResult())

            if self.analyticalField.type == "scalarExpression" and not currentResults.shape[2] == 1:
                raise Exception("Cannot map scalar value to {}-dimensional result.".format(currentResults.shape[2]))

            elementList = self.fieldOutput.elSet

            for i1, element in enumerate(elementList):
                coordinatesAtCenter = element.getCoordinatesAtCenter()
                # for i2, quadraturePoint in enumerate(self.fieldOutput.quadraturePoints):
                #    currentResults[i1][i2] = self.analyticalField.evaluateAtCoordinates(coordinatesAtCenter)
                currentResults[i1] = self.analyticalField.evaluateAtCoordinates(coordinatesAtCenter)

            self.fieldOutput.setResults(currentResults)
