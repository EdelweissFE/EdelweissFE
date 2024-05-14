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

import textwrap
from collections import defaultdict
from warnings import warn

from edelweissfe.config.stepactions import stepActionFactory
from edelweissfe.journal.journal import Journal
from edelweissfe.models.femodel import FEModel
from edelweissfe.steps.adaptivestep import AdaptiveStep
from edelweissfe.steps.base.stepbase import StepBase
from edelweissfe.utils.fieldoutput import FieldOutputController


class StepActionDefinition:
    def __init__(self, name: str, module: str, kwargs: dict):
        self.name = name
        self.module = module
        self.kwargs = kwargs


class StepDefinition:
    def __init__(
        self,
        stepType: str,
        stepOptions: dict,
        stepActionDefinitions: list[StepActionDefinition],
    ):
        self.type = stepType
        self.stepOptions = stepOptions
        self.stepActionDefinitions = stepActionDefinitions


class StepManager:
    """This manager for convenience parses all step defintions for the simulation
    and calls the respective modules, which generate (or update) StepActions based on
    computed results, model info and job information.

    Parameters
    ----------
    inputfile
        The inputfile containing the step defintions.
    """

    identification = "StepManager"

    def __init__(
        self,
    ):
        self.stepActions = defaultdict(dict)
        self.stepDefinitions = []

    def enqueueStepDefinition(self, stepDefinition: StepDefinition):
        """Enqueue a step definition.

        Parameters
        ----------
        stepDefinition
            The StepDefinition containing the step type, step options and list of StepActionDefinition.
        """
        self.stepDefinitions.append(stepDefinition)

    def dequeueStep(
        self,
        jobInfo: dict,
        model: FEModel,
        fieldOutputController: FieldOutputController,
        journal: Journal,
        solvers: dict,
        outputManagers: list,
    ) -> StepBase:
        """Dequeue the next step.

        Parameters
        ----------
        jobInfo
            A dictionary containing information on the job.
        model
            The model tree.
        stepActions
            A dictionary containing already existing step actions.
        fieldOutputController
            The field output controller.
        journal
            The journal instance for logging.
        outputManagers
            The OutputManagers used.
        stepActions
            The collection of actions for this step.

        Returns
        -------
        StepActionBase
            The next Step.
        """

        def printActionDefinition(intro, options):
            for line in textwrap.wrap(
                intro + " [" + ", ".join(("{:}={:}".format(k, v) for k, v in options.items())) + "]",
                subsequent_indent=" " * (len(intro) + 1),
            ):
                journal.message(
                    line,
                    self.identification,
                    2,
                )

        for stepNumber, stepDefinition in enumerate(self.stepDefinitions):
            actionDefinitionsInThisStep = []

            journal.message(
                "StepAction definitions:",
                self.identification,
                1,
            )

            for action in stepDefinition.stepActionDefinitions:
                if action.name in actionDefinitionsInThisStep:
                    raise Exception(
                        "Warning: StepAction {:} has multiple definitions in step {:}".format(action.name, stepNumber)
                    )
                actionDefinitionsInThisStep.append(action.name)

                if action.name in self.stepActions[action.module]:
                    self.stepActions[action.module][action.name].updateStepAction(
                        action.kwargs, jobInfo, model, fieldOutputController, journal
                    )
                    printActionDefinition('Updating "{:}"'.format(action.name), action.kwargs)

                else:
                    printActionDefinition('Creating "{:}"'.format(action.name), action.kwargs)

                    self.stepActions[action.module][action.name] = stepActionFactory(action.module)(
                        action.name,
                        action.kwargs,
                        jobInfo,
                        model,
                        fieldOutputController,
                        journal,
                    )

            try:
                solverName = stepDefinition.stepOptions.pop("solver")
            except KeyError:
                # raise KeyError("Step definition missing solver option.")
                warn(
                    "Step definition missing solver option.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                solverName = "default"

            try:
                solver = solvers[solverName]
            except KeyError:
                mssg = f"No definition found for solver {solverName}."
                availableSolvers = [key for key in solvers.keys() if not key == "default"]
                if availableSolvers:
                    mssg += " Available solvers: " + ", ".join(availableSolvers)
                else:
                    mssg += " Define solver using *solver keyword."
                raise KeyError(mssg)

            yield AdaptiveStep(
                stepNumber,
                model,
                fieldOutputController,
                journal,
                jobInfo,
                solver,
                outputManagers,
                self.stepActions,
                **stepDefinition.stepOptions,
            )
