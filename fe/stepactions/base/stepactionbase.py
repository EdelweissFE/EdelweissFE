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
"""
Created on Sat May  6 20:05:41 2017

@author: Matthias Neuner
"""

from abc import ABC, abstractmethod
from fe.utils.fieldoutput import FieldOutputController
from fe.journal.journal import Journal
from fe.utils.dofmanager import DofVector


class StepActionBase(ABC):
    """This is the abstract base class for all step actions managers.
    User defined step actions must implement the abstract methods."""

    identification = "StepActionBase"

    @abstractmethod
    def __init__(
        self,
        name: str,
        definition: dict,
        jobInfo: dict,
        modelInfo: dict,
        fieldOutputController: FieldOutputController,
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
        modelInfo
            A dictionary containing the model tree.
        fieldOutputController
            The fieldput controlling object.
        journal
            The journal object for logging.
        """

        pass

    @abstractmethod
    def updateStepAction(
        self,
        name: str,
        definition: dict,
        jobInfo: dict,
        modelInfo: dict,
        fieldOutputController: FieldOutputController,
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
        modelInfo
            A dictionary containing the model tree.
        fieldOutputController
            The fieldput controlling object.
        journal
            The journal object for logging.
        """

        pass

    @abstractmethod
    def finishStep(self, U: DofVector, P: DofVector):
        """Is called when a step successfully finished.

        Parameters
        ----------
        U
            The current solution vector.
        P
            The current reaction force vector.
        """
        pass
