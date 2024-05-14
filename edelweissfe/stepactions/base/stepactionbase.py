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

from abc import ABC

from edelweissfe.journal.journal import Journal
from edelweissfe.models.femodel import FEModel
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.utils.fieldoutput import FieldOutputController


class StepActionBase(ABC):
    """This is the abase class for all step actions.
    User defined step actions can override the methods.

    Parameters
    ----------
    name
        The name of this step action.
    definition
        A dictionary containing the options for this step action.
    jobInfo
        A dictionary containing the information about the job.
    model
        A dictionary containing the model tree.
    fieldOutputController
        The fieldput controlling object.
    journal
        The journal object for logging.
    """

    def __init__(
        self,
        name: str,
        definition: dict,
        jobInfo: dict,
        model: FEModel,
        fieldOutputController: FieldOutputController,
        dofmanager,
        journal: Journal,
    ):
        pass

    def updateStepAction(
        self,
        definition: dict,
        jobInfo: dict,
        model: FEModel,
        fieldOutputController: FieldOutputController,
        dofmanager,
        journal: Journal,
    ):
        """Is called when an updated definition is present for a new step.

        Parameters
        ----------
        name
            The name of this step action.
        definition
            A dictionary containing the options for this step action.
        jobInfo
            A dictionary containing the information about the job.
        model
            A dictionary containing the model tree.
        fieldOutputController
            The fieldput controlling object.
        journal
            The journal object for logging.
        """

    def applyAtStepStart(self, model):
        """Is called when a step starts.

        Parameters
        ----------
        U
            The current solution vector.
        P
            The current reaction force vector.
        """

    def applyAtStepEnd(self, model):
        """Is called when a step successfully finished.

        Parameters
        ----------
        U
            The current solution vector.
        P
            The current reaction force vector.
        """

    def applyAtIncrementStart(self, model, timeStep: TimeStep):
        """Is called when a step increment starts.

        Parameters
        ----------
        U_n
            The current converged solution vector at start of the increment.
        P
            The current reaction force vector.
        increment
            The defintion of the time increment.
        """
