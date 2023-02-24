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
# Created on Tue Jan  17 19:10:42 2017

# @author: Matthias Neuner
"""This is the main module of EdelweissFE.

Heart is the ``*job`` keyword, which defines the spatial dimension
A ``*job`` definition consists of multiple ``*steps``, associated with that job.
"""

from collections import OrderedDict, defaultdict
from fe.config import analyticalfields
from fe.config.phenomena import domainMapping
from fe.config.generators import getGeneratorFunction
from fe.config.stepactions import stepActionFactory
from fe.config.outputmanagers import getOutputManagerClass
from fe.config.solvers import getSolverByName
from fe.utils.fieldoutput import FieldOutputController
from fe.utils.misc import convertAssignmentsToStringDictionary, splitLineAtCommas
from fe.utils.plotter import Plotter
from fe.utils.exceptions import StepFailed
from fe.utils.dofmanager import DofManager, DofVector
from fe.variables.scalarvariable import ScalarVariable
from fe.config.configurator import loadConfiguration, updateConfiguration
from fe.journal.journal import Journal
from fe.utils.caseinsensitivedict import CaseInsensitiveDict
from fe.utils.abqmodelconstructor import AbqModelConstructor
from time import time as getCurrentTime

from fe.stepactions.base.stepactionbase import StepActionBase


def gatherStepActions(
    step: dict,
    jobInfo: dict,
    model: dict,
    time: float,
    U: DofVector,
    P: DofVector,
    stepActions: dict[str, StepActionBase],
    fieldOutputController: FieldOutputController,
    journal: Journal,
) -> dict[str, StepActionBase]:
    """Parses all the defined actions for the current step,
    and calls the respective modules, which generate step-actions based on
    computed results, model info and job information.
    The stepi actions are stored in a dictionary, which is handed to
    solveStep() in the feCore main routine afterwards.
    The step action modules decide if old stepaction definitions are
    overwritten or extended. Returns a dictionary with keys as defined in
    stepactions.

    Parameters
    ----------
    step
        The step active for which the actions are gathered.
    jobInfo
        A dictionary containing information on the job.
    model
        A dictionary containing the model tree.
    time
        The current time of the simulation.
    U
        The current solution vector.
    P
        The current reaction vector.
    stepActions
        A dictionary containing already existing step actions.
    fieldOutputController
        The field output controller.
    journal
        The journal instance for logging.

    Returns
    -------
    dict[str,StepActionBase]
        The updated dictionary of step actions.
    """

    for dataline in step["data"]:
        actionType, *definition = splitLineAtCommas(dataline)
        options = convertAssignmentsToStringDictionary(definition)

        module = actionType.lower()
        moduleName = options.get("name", options.get("category", module))

        if moduleName in stepActions[module]:
            stepActions[module][moduleName].updateStepAction(
                moduleName, options, jobInfo, model, fieldOutputController, journal
            )
            journal.message('Stepaction "{:}" will be updated'.format(moduleName), "stepActionManager", 1)
        else:
            stepActions[module][moduleName] = stepActionFactory(module)(
                moduleName, options, jobInfo, model, fieldOutputController, journal
            )

    return stepActions


def finiteElementSimulation(
    inputfile: dict, verbose: bool = False, suppressPlots: bool = False
) -> tuple[bool, DofVector, DofVector, FieldOutputController]:
    """This is core function of the finite element analysis.
    Based on the keyword ``*job``, the finite element model is defined.

    It assembles
     * the information on the job
     * the model tree
     * steps
     * field outputs
     * output managers

    and controls the respective solver based on the defined simulation steps.
    For each step, the step-actions (dirichlet, nodeforces) are collected by
    external modules.

    Parameters
    ----------
    inputfile
        The input file in dictionary form.
    verbose
        Be verbose during the simulation.
    suppressPlots
        Suppress plots at the end of simulation for batch runs.

    Returns
    -------
    tuple
        A tuple containing
            - Truth value of success
            - The solution vector.
            - The reaction vector.
            - The fieldoutput controller containing all results.

    """

    identification = "feCore"

    journal = Journal(verbose=verbose)

    # create job dictionary from input file like e.g.
    job = inputfile["*job"][0]

    modelDomain = job["domain"]
    domainSize = domainMapping[modelDomain]

    model = {
        "nodes": {},
        "elements": {},
        "nodeSets": {},
        "elementSets": {},
        "surfaces": {},
        "constraints": {},
        "materials": {},
        "analyticalFields": {},
        "scalarVariables": [],
        "domainSize": domainSize,
    }

    # call individual optional model generators
    for generatorDefinition in inputfile["*modelGenerator"]:
        if generatorDefinition.get("executeAfterManualGeneration", False):
            continue
        gen = generatorDefinition["generator"]
        model = getGeneratorFunction(gen)(generatorDefinition, model, journal)

    # the standard 'Abaqus like' model generator is invoked unconditionally, and it has direct access to the inputfile
    abqModelConstructor = AbqModelConstructor(journal)
    model = abqModelConstructor.createGeometryFromInputFile(model, inputfile)
    model = abqModelConstructor.createMaterialsFromInputFile(model, inputfile)
    model = abqModelConstructor.createConstraintsFromInputFile(model, inputfile)
    model = abqModelConstructor.createAnalyticalFieldsFromInputFile(model, inputfile)
    model = abqModelConstructor.createSectionsFromInputFile(model, inputfile)

    # call individual optional model generators,
    for generatorDefinition in inputfile["*modelGenerator"]:
        if not generatorDefinition.get("executeAfterManualGeneration", False):
            continue
        gen = generatorDefinition["generator"]
        model = getGeneratorFunction(gen)(generatorDefinition, model, journal)

    # we may have additional scalar degrees of freedom, not associated with any node (e.g, lagrangian multipliers of constraints)
    for constraintName, constraint in model["constraints"].items():
        nAdditionalScalarVariables = constraint.getNumberOfAdditionalNeededScalarVariables()
        journal.message(
            "Constraint {:} requests {:} additional scalar dofs".format(constraintName, nAdditionalScalarVariables),
            identification,
            0,
        )

        scalarVariables = [ScalarVariable() for i in range(nAdditionalScalarVariables)]
        model["scalarVariables"] += scalarVariables

        constraint.assignAdditionalScalarVariables(scalarVariables)

    # create total number of dofs and orderedDict of fieldType and associated numbered dofs
    dofManager = DofManager(model)

    journal.message("total size of eq. system: {:}".format(dofManager.nDof), identification, 0)
    journal.printSeperationLine()

    jobName = job.get("name", "")
    time = job.get("startTime", 0.0)
    # store job info
    jobInfo = {"dofManager": dofManager, "computationTime": 0.0}

    # add or update additional job info such as inputfile, domain, name, tolerances
    jobInfo.update(job)
    jobInfo = loadConfiguration(jobInfo)
    for updateConfig in inputfile["*updateConfiguration"]:
        updateConfiguration(updateConfig, jobInfo, journal)

    plotter = Plotter(journal, inputfile)

    # generate an instance of the desired solver
    Solver = getSolverByName(job.get("solver", "NIST"))
    solver = Solver(jobInfo, journal)
    U, P = solver.initialize()
    model["solver"] = solver

    fieldOutputController = FieldOutputController(model, inputfile, journal)

    # collect all output managers in a list of objects
    outputmanagers = []
    for outputDef in inputfile["*output"]:
        OutputManager = getOutputManagerClass(outputDef["type"].lower())
        managerName = outputDef.get("name", jobName + outputDef["type"])
        definitionLines = outputDef["data"]
        outputmanagers.append(
            OutputManager(managerName, definitionLines, jobInfo, model, fieldOutputController, journal, plotter)
        )

    stepActions = defaultdict(CaseInsensitiveDict)

    fieldOutputController.initializeJob(time, U, P)

    for man in outputmanagers:
        man.initializeSimulation(model)
    try:
        for stepNumber, step in enumerate(inputfile["*step"]):
            try:
                stepActions = gatherStepActions(
                    step, jobInfo, model, time, U, P, stepActions, fieldOutputController, journal
                )

                for modelUpdate in stepActions["modelupdate"].values():
                    model = modelUpdate.updateModel(model, fieldOutputController, journal)

                fieldOutputController.initializeStep(step, stepActions)
                for manager in outputmanagers:
                    manager.initializeStep(step, stepActions)

                # solve the step
                tic = getCurrentTime()
                success, U, P, time = solver.solveStep(
                    stepNumber, step, time, stepActions, model, U, P, fieldOutputController, outputmanagers
                )
                toc = getCurrentTime()

                stepTime = toc - tic
                jobInfo["computationTime"] += stepTime

                if not success:
                    raise StepFailed()

                journal.printTable(
                    [
                        ("Step computation time", "{:10.4f}s".format(stepTime)),
                    ],
                    identification,
                    level=0,
                )

            finally:
                fieldOutputController.finalizeStep(U, P)
                for manager in outputmanagers:
                    manager.finalizeStep(U, P, time)

    except KeyboardInterrupt:
        print("")
        journal.errorMessage("Interrupted by user", identification)
        success = True

    except StepFailed:
        success = False
        journal.errorMessage("Step not finished", identification)

    finally:
        journal.printTable(
            [
                ("Job computation time", "{:10.4f}s".format(jobInfo["computationTime"])),
            ],
            identification,
            level=0,
            printHeaderRow=False,
        )
        # let all output managers finalize the job
        fieldOutputController.finalizeJob(U, P)
        for manager in outputmanagers:
            manager.finalizeJob(
                U,
                P,
            )

        plotter.finalize()
        if not suppressPlots:
            plotter.show()

    return success, U, P, fieldOutputController
