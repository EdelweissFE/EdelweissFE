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

from edelweissfe.journal.journal import Journal
from edelweissfe.models.femodel import FEModel
from edelweissfe.steps.base.stepbase import StepBase
from edelweissfe.timesteppers.adaptivetimestepper import AdaptiveTimeStepper
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.utils.caseinsensitivedict import CaseInsensitiveDict
from edelweissfe.utils.fieldoutput import FieldOutputController


class AdaptiveStep(StepBase):
    """
    A standard adaptive incremental step to be used in nonlinear simulations.

    Parameters
    ----------
    number
        The number of this step. For information purposes only.
    starTime
        The start time of the step.
    definition
        A dictionary holding key/value pairs for defintion
    stepActions
        The collection of actions for this step.
    journal
        The journal object for logging.
    """

    defaultStartInc = 1.0
    defaultMaxInc = 1.0
    defaultMinInc = 1e-4
    defaultMaxNumInc = 1000
    defaultMaxIter = 10
    defaultCriticalIter = 5
    defaultMaxGrowingIter = 10

    def __init__(
        self,
        number: int,
        model: FEModel,
        fieldOutputController: FieldOutputController,
        journal: Journal,
        jobInfo: dict,
        solver,
        outputManagers: list,
        stepActions: dict,
        **kwargs
    ):
        kwargs = CaseInsensitiveDict(kwargs)
        self.number = number  #: The (unique) number of the step.

        self.model = model
        self.fieldOutputController = fieldOutputController
        self.journal = journal
        self.solver = solver
        self.outputManagers = outputManagers

        self.length = kwargs.get("stepLength", 1.0)  #: The durcation of the step.
        self.startIncrementSize = kwargs.get(
            "startInc", self.defaultStartInc
        )  #: The initial fraction of the step to be computed.
        self.maxIncrementSize = kwargs.get(
            "maxInc", self.defaultMaxInc
        )  #: The maximal fraction of the step to be computed.
        self.minIncrementSize = kwargs.get(
            "minInc", self.defaultMinInc
        )  #: The minimal fraction of the step to be computed.
        self.maxNumberIncrements = kwargs.get(
            "maxNumInc", self.defaultMaxNumInc
        )  #: The maximal number of increments allowed.
        self.maxIter = kwargs.get("maxIter", self.defaultMaxIter)  #: The maximal number of iterations allowed.
        self.criticalIter = kwargs.get(
            "criticalIter", self.defaultCriticalIter
        )  #: The number of critical iterations after which the next increment is reduced.
        self.maxGrowIter = kwargs.get(
            "maxGrowIter", self.defaultMaxGrowingIter
        )  #: The number of residual growths before the increment is discarded.

        self.incrementGenerator = AdaptiveTimeStepper(
            model.time,
            self.length,
            self.startIncrementSize,
            self.maxIncrementSize,
            self.minIncrementSize,
            self.maxNumberIncrements,
            journal,
        )

        self.actions = stepActions

    def solve(
        self,
    ):
        model = self.model
        fieldOutputController = self.fieldOutputController
        journal = self.journal
        outputManagers = self.outputManagers

        try:
            for modelUpdate in self.actions["modelupdate"].values():
                model = modelUpdate.updateModel(model, fieldOutputController, journal)

            fieldOutputController.initializeStep(self)
            for manager in outputManagers:
                manager.initializeStep(self)

            self.solver.solveStep(self, model, fieldOutputController, outputManagers)

        finally:
            fieldOutputController.finalizeStep()
            for manager in outputManagers:
                manager.finalizeStep()

    def getTimeStep(
        self,
    ) -> TimeStep:
        return self.incrementGenerator.generateTimeStep()

    def discardAndChangeIncrement(self, cutbackFactor: float):
        return self.incrementGenerator.discardAndChangeIncrement(cutbackFactor)

    def preventIncrementIncrease(
        self,
    ):
        return self.incrementGenerator.preventIncrementIncrease()
