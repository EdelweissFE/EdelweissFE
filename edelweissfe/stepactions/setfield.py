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

import numpy as np

from edelweissfe.stepactions.base.stepactionbase import StepActionBase

"""
Set a field (via fieldOutput) to a predefined value.
"""

documentation = {
    "fieldOutput": "Field output to be set",
    "type": "'Const' or 'analyticalField'",
    "value": "Scalar value if type 'const'; name of analyticalField if type 'analyticalField'",
}


class StepAction(StepActionBase):
    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name
        self.active = True
        self.journal = journal
        self.fieldOutputName = action["fieldOutput"]
        self.fieldOutput = fieldOutputController.fieldOutputs[self.fieldOutputName]
        self.type = action["type"]
        self.value = action["value"]

        if self.type not in ["uniform", "analyticalField"]:
            raise Exception("Invalid type: {}".format(self.type))

        if self.type == "analyticalField":
            self.analyticalField = model.analyticalFields[self.value]

    def applyAtStepEnd(self, model, stepMagnitude=None):
        self.active = False

    def updateStepAction(self, action, jobInfo, model, fieldOutputController, journal):
        self.active = True

    def applyAtStepStart(
        self,
        model,
    ):
        if not self.active:
            return

        if self.type == "uniform":
            currentResults = np.zeros_like(self.fieldOutput.getLastResult())
            newResult = np.fromstring(self.value, float, sep=",")

            if currentResults.ndim == 2:
                currentResults = np.expand_dims(currentResults, 0)

            if not currentResults.shape[-1] == newResult.shape[-1]:
                self.journal.errorMessage(
                    f"Dimension mismatch. Result '{self.fieldOutputName}' has length {currentResults.shape[-1]} but value has length {newResult.shape[-1]}",  # noqa: E501
                    self.name,
                )
                raise Exception

            currentResults[:] = newResult
            self.journal.message(
                "Setting field {:} to uniform value {:}".format(self.fieldOutputName, self.value),
                self.name,
            )
            self.fieldOutput.setResults(currentResults)

        if self.type == "analyticalField":
            currentResults = np.zeros_like(self.fieldOutput.getLastResult())

            if self.analyticalField.type == "scalarExpression" and not currentResults.shape[2] == 1:
                self.journal.errorMessage(f"Cannot map scalar value to {currentResults.shape[2]}-dimensional result.")
                raise Exception

            elementList = self.fieldOutput.associatedSet

            for i1, element in enumerate(elementList):
                coordinatesAtQuadraturePoints = element.getCoordinatesAtQuadraturePoints()

                currentResults[i1] = self.analyticalField.evaluateAtCoordinates(coordinatesAtQuadraturePoints)

            self.fieldOutput.setResults(currentResults)
