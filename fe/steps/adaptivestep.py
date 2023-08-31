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

from collections import defaultdict
from fe.steps.base.stepbase import StepBase
from fe.utils.incrementgenerator import IncrementGenerator
from fe.journal.journal import Journal


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
        self, number: int, startTime: float, definition: defaultdict, stepActions: dict, jobInfo: dict, journal: Journal
    ):
        self.number = number
        self.length = definition.get("stepLength", 1.0)
        self.startIncrementSize = definition.get("startInc", self.defaultStartInc)
        self.maxIncrementSize = definition.get("maxInc", self.defaultMaxInc)
        self.minIncrementSize = definition.get("minInc", self.defaultMinInc)
        self.maxNumberIncrements = definition.get("maxNumInc", self.defaultMaxNumInc)
        self.maxIter = definition.get("maxIter", self.defaultMaxIter)
        self.criticalIter = definition.get("criticalIter", self.defaultCriticalIter)
        self.maxGrowIter = definition.get("maxGrowIter", self.defaultMaxGrowingIter)

        self.incrementGenerator = IncrementGenerator(
            startTime,
            self.length,
            self.startIncrementSize,
            self.maxIncrementSize,
            self.minIncrementSize,
            self.maxNumberIncrements,
            journal,
        )

        self.actions = stepActions

        self.solverName = definition.get("solver", None)

    def solve(self, solvers, model, fieldOutputController, outputManagers, journal):
        if self.solverName in solvers:
            solver = solvers[self.solverName]
        else:
            from warnings import warn

            warn(
                "Warning, not defining a solver for a step is deprecated; assign a solver using the solver= option",
                DeprecationWarning,
                stacklevel=2,
            )
            solver = solvers["default"]

        success = False

        try:
            for modelUpdate in self.actions["modelupdate"].values():
                model = modelUpdate.updateModel(model, fieldOutputController, journal)

            fieldOutputController.initializeStep(self)
            for manager in outputManagers:
                manager.initializeStep(self)

            success, model = solver.solveStep(self, model, fieldOutputController, outputManagers)

        finally:
            fieldOutputController.finalizeStep(model)
            for manager in outputManagers:
                manager.finalizeStep(model)

        return success, model

    def getTimeIncrement(
        self,
    ):
        return self.incrementGenerator.generateIncrement()

    def discardAndChangeIncrement(self, cutbackFactor: float):
        return self.incrementGenerator.discardAndChangeIncrement(cutbackFactor)

    def preventIncrementIncrease(
        self,
    ):
        return self.incrementGenerator.preventIncrementIncrease()
