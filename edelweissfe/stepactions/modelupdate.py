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

from edelweissfe.stepactions.base.stepactionbase import StepActionBase
from edelweissfe.utils.math import execModelAccessibleExpression

"""This step action may be used for updating something in the model at the beginning
of a step.
"""

documentation = {"update": "Model accessible, executable expression"}


class StepAction(StepActionBase):
    def __init__(self, name, options, jobInfo, model, fieldOutputController, journal):
        self.name = name
        self.updateStepAction(options, jobInfo, model, fieldOutputController, journal)

    def applyAtStepEnd(self, model):
        """By default, this action is only executed once."""

        self.active = False

    def updateStepAction(self, options, jobInfo, model, fieldOutputController, journal):
        """Update the expression, and set the action active again."""

        self.updateExpression = options["update"]
        self.active = True

    def updateModel(self, model, fieldOutputController, journal):
        """Update the model based on an executable provided Python expression."""

        if not self.active:
            return

        journal.message("Updating model: {:}".format(self.updateExpression), self.name)
        execModelAccessibleExpression(
            self.updateExpression,
            model,
            fieldOutputs=fieldOutputController.fieldOutputs,
        )
        return model
