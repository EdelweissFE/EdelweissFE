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

from abc import ABC, abstractmethod

from edelweissfe.journal.journal import Journal
from edelweissfe.models.femodel import FEModel
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.utils.fieldoutput import FieldOutputController
from edelweissfe.utils.plotter import Plotter


class OutputManagerBase(ABC):
    """This is the abstract base class for all output managers.
    User defined output managers must implement the abstract methods.

    Parameters
    ----------
    name
        The name of this output manager.
    definitionLines
        The dictionary containing the definition of the output manager.
    model
        A dictionary containing the model tree.
    fieldOutputController
        The field output contoller instance.
    journal
        The journal instance for logging.
    plotter
        The plotter instance for plotting.
    """

    identification = "OutputManagerBase"

    @abstractmethod
    def __init__(
        self,
        name: str,
        definitionLines: dict,
        model: FEModel,
        fieldOutputController: FieldOutputController,
        journal: Journal,
        plotter: Plotter,
    ):
        pass

    @abstractmethod
    def updateDefinition(self, **kwargs: dict):
        pass

    @abstractmethod
    def initializeJob(self):
        """Initalize the output manager at the beginning of a step.

        Parameters
        ----------
        """

    @abstractmethod
    def initializeStep(self, step: dict):
        """Initalize the output manager at the beginning of a step.

        Parameters
        ----------
        step
            A dictionary containing the step definition.
        """

    @abstractmethod
    def finalizeIncrement(self, timeStep: TimeStep, **kwargs):
        """Finalize the output at the end of a time increment.

        Parameters
        ----------
        U
            The initial solution vector.
        P
            The initial reaction vector.
        timeStep
            The time step.
        **kwargs
            Keyword arguments.
        """

    @abstractmethod
    def finalizeFailedIncrement(self, **kwargs):
        """Finalize the output at the end of a time increment.

        Parameters
        ----------
        **kwargs
            Keyword arguments.
        """

    @abstractmethod
    def finalizeStep(
        self,
    ):
        """Finalize the output the end of a step."""

    @abstractmethod
    def finalizeJob(
        self,
    ):
        """Finalize the output at the end of a job.

        Parameters
        ----------
        U
            The final solution vector.
        P
            The final reaction vector.
        """
