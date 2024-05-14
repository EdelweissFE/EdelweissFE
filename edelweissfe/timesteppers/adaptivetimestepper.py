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
# Created on Sat Jan  21 12:18:10 2017

from edelweissfe.journal.journal import Journal
from edelweissfe.timesteppers.timestep import TimeStep

# @author: Matthias
from edelweissfe.utils.exceptions import ReachedMaxIncrements, ReachedMinIncrementSize


class AdaptiveTimeStepper:
    identification = "AdaptiveTimeStepper"

    def __init__(
        self,
        currentTime: float,
        stepLength: float,
        startIncrement: float,
        maxIncrement: float,
        minIncrement: float,
        maxNumberIncrements: int,
        journal: Journal,
    ):
        """
        An increment generator for incremental-iterative simulations.

        Implementation as generator class.

        Parameters
        ----------
        currentTime
            The current (start) time.
        stepLength
            The total length of the step.
        startIncrement
            The size of the start increment.
        maxIncrement
            The maximum size of an increment.
        minIncrement
            The minimum size of an increment.
        maxNumberIncrements
            The maximum number of allowed increments.
        journal
            The journal instance for logging purposes.
        """

        self.nPassedGoodIncrements = int(0)
        self.totalIncrements = int(0)
        self.startIncrement = startIncrement
        self.maxIncrement = maxIncrement
        self.minIncrement = minIncrement
        self.maxNumberIncrements = maxNumberIncrements

        self.finishedStepProgress = 0.0
        self.increment = min(startIncrement, maxIncrement)
        self.allowedToIncreasedNext = True

        self.currentTime = currentTime
        self.stepLength = stepLength
        self.dT = 0.0
        self.journal = journal

    def generateTimeStep(self) -> TimeStep:
        """
        Generate the next increment.

        Returns
        -------
        TimeStep
            The current time step.
        """

        # zero increment; return value for first function call
        yield TimeStep(0, 0.0, 0.0, 0.0, 0.0, self.currentTime)

        while self.finishedStepProgress < (1.0 - 1e-15):
            if self.totalIncrements >= self.maxNumberIncrements:
                self.journal.errorMessage("Reached maximum number of increments", self.identification)
                raise ReachedMaxIncrements()

            if (self.nPassedGoodIncrements >= 3) and self.allowedToIncreasedNext:
                self.increment *= 1.1
                if self.increment > self.maxIncrement:
                    self.increment = self.maxIncrement
            self.allowedToIncreasedNext = True

            remainder = 1.0 - self.finishedStepProgress
            if remainder < self.increment:
                self.increment = remainder

            dT = self.stepLength * self.increment
            self.finishedStepProgress += self.increment
            endTimeOfIncrementInStep = self.stepLength * self.finishedStepProgress
            endTimeOfIncrementInTotal = self.currentTime + endTimeOfIncrementInStep

            self.totalIncrements += 1
            self.nPassedGoodIncrements += 1

            yield TimeStep(
                self.totalIncrements,
                self.increment,
                self.finishedStepProgress,
                dT,
                endTimeOfIncrementInStep,
                endTimeOfIncrementInTotal,
            )

    def preventIncrementIncrease(
        self,
    ):
        """May be called before an increment is requested, to prevent from
        automatically increasing, e.g. in case of bad convergency."""

        self.allowedToIncreasedNext = False

    def discardAndChangeIncrement(self, scaleFactor: float):
        """Change increment size between minIncrement and
        maxIncrement by a given scale factor.

        Parameters
        ----------
        scaleFactor
            The factor for scaling based on the previous increment.
        """

        if self.increment == self.minIncrement:
            self.journal.errorMessage("Cannot reduce increment size", self.identification)
            raise ReachedMinIncrementSize()

        if self.finishedStepProgress == 0.0:
            self.journal.errorMessage("Failed zero increment", self.identification)
            raise ReachedMinIncrementSize()

        self.finishedStepProgress -= self.increment
        newIncrement = self.increment * scaleFactor
        if newIncrement > self.maxIncrement:
            self.increment = self.maxIncrement
        elif newIncrement < self.minIncrement:
            self.increment = self.minIncrement
        else:
            self.increment = newIncrement

        self.journal.message(
            "Cutback to increment size {:}".format(self.increment),
            self.identification,
            2,
        )
        self.totalIncrements -= 1
        self.nPassedGoodIncrements = 0
