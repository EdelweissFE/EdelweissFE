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

from fe.stepactions.base.stepactionbase import StepActionBase
from fe.utils.fieldoutput import FieldOutputController
from fe.journal.journal import Journal
from fe.utils.misc import convertAssignmentsToStringDictionary, splitLineAtCommas
from fe.config.stepactions import stepActionFactory
from fe.models.femodel import FEModel
from collections import OrderedDict, defaultdict
from fe.utils.caseinsensitivedict import CaseInsensitiveDict
from fe.steps.adaptivestep import AdaptiveStep
from fe.steps.base.stepbase import StepBase
import textwrap


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

    def __init__(self, inputfile):
        self.stepActions = defaultdict(CaseInsensitiveDict)
        self.inputfile = inputfile
        self.steps = 0

    def getStep(
        self,
        jobInfo: dict,
        model: FEModel,
        fieldOutputController: FieldOutputController,
        journal: Journal,
    ) -> StepBase:
        """Get the next step.

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

        Returns
        -------
        StepActionBase
            The next Step.
        """

        def printActionDefinition(intro, options):
            for l in textwrap.wrap(
                intro + " [" + ", ".join(("{:}={:}".format(k, v) for k, v in options.items())) + "]",
                subsequent_indent=" " * (len(intro) + 1),
            ):
                journal.message(
                    l,
                    self.identification,
                    2,
                )

        for stepNumber, stepDefinition in enumerate(self.inputfile["*step"]):
            actionDefinitionsInThisStep = []

            journal.message(
                "StepAction definitions:",
                self.identification,
                1,
            )

            for dataline in stepDefinition["data"]:
                actionType, *definition = splitLineAtCommas(dataline)
                options = convertAssignmentsToStringDictionary(definition)

                module = actionType.lower()
                moduleName = options.get("name", options.get("category", module))

                if moduleName in actionDefinitionsInThisStep:
                    raise Exception(
                        "Warning: StepAction {:} has multiple definitions in step {:}".format(moduleName, stepNumber)
                    )

                actionDefinitionsInThisStep.append(moduleName)

                stepActions = self.stepActions

                if moduleName in stepActions[module]:
                    stepActions[module][moduleName].updateStepAction(
                        moduleName, options, jobInfo, model, fieldOutputController, journal
                    )
                    printActionDefinition('Updating "{:}"'.format(moduleName), options)

                else:
                    printActionDefinition('Creating "{:}"'.format(moduleName), options)

                    stepActions[module][moduleName] = stepActionFactory(module)(
                        moduleName, options, jobInfo, model, fieldOutputController, journal
                    )

            yield AdaptiveStep(stepNumber, model.time, stepDefinition, stepActions, jobInfo, journal)
