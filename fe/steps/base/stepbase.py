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

from fe.journal.journal import Journal
from fe.models.femodel import FEModel
from fe.utils.fieldoutput import FieldOutputController
from collections import defaultdict


class StepBase:
    """
    This is a simulation step.

    It has a specific runtime, and it holds StepActions to executed.


    Parameters
    ----------
    number
        The number of this step. For information purposes only.
    startTime
        The start time of the step.
    definition
        A dictionary holding key/value pairs for defintion
    stepActions
        The collection of actions for this step.
    journal
        The journal object for logging.
    """

    def __init__(self, number: int, startTime: float, definition: defaultdict, stepActions: dict, journal):
        pass

    def solve(
        self, solvers: dict, model: FEModel, fieldOutputController: FieldOutputController, outputManagers: dict
    ) -> tuple[bool, FEModel]:
        """
        Let a step be solved.

        Parameters
        ----------
        solvers
            The instances of solvers available to this step.
        model
            The model to be solved.
        fieldOutputController
            The FieldOutputController instance for post processing of results.
        outputManagers
            The OutputManagers used.
        journal
            The journal instance for logging.

        Returns
        -------
        tuple[bool, FEModel]
            The tuple containing the truth value of successs and the updated model.
        """

        pass
