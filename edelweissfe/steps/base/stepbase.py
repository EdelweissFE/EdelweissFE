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
"""One or more Steps form the simulation. Steps are executed in consecutive order,
have a distinct (physical) runtime, and may contain multiple StepActions.
Subsequent Steps inherit StepActions, and they may be updated."""

from edelweissfe.journal.journal import Journal
from edelweissfe.models.femodel import FEModel
from edelweissfe.utils.fieldoutput import FieldOutputController


class StepBase:
    """
    This is a simulation step.

    It has a specific runtime, and it holds StepActions to executed.

    Parameters
    ----------
    number
        The number of this step. For information purposes only.
    model
        The current state of the model.
    fieldOutputController
        The FieldOutputController instance for processing results.
    journal
        The Journal instance for logging purposes.
    jobInfo
        Additional information about the job
    solvers
        The instances of solvers available to this step.
    outputManagers
        The OutputManagers used.
    stepActions
        The collection of actions for this step.
    **kwargs
        Additional options for the step.
    """

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
        **kwargs,
    ):
        pass

    def solve(
        self,
    ) -> FEModel:
        """
        Let a step be solved.

        Parameters
        ----------

        Returns
        -------
        FEModel
            The updated model.
        """
