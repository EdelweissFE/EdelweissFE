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

from time import time as getCurrentTime

from edelweissfe.config.configurator import loadConfiguration, updateConfiguration
from edelweissfe.config.phenomena import domainMapping
from edelweissfe.config.solvers import getSolverByName
from edelweissfe.helpers.inputfilehelpers import (
    createFieldOutputFromInputFile,
    createOutputManagersFromInputFile,
    createPlotterFromInputFile,
    createSolversFromInputFile,
    createStepManagerFromInputFile,
    fillFEModelFromInputFile,
)
from edelweissfe.journal.journal import Journal
from edelweissfe.models.femodel import FEModel, printPrettyModelSummary
from edelweissfe.utils.exceptions import StepFailed
from edelweissfe.utils.fieldoutput import FieldOutputController


def finiteElementSimulation(
    inputfile: dict, verbose: bool = False, suppressPlots: bool = False
) -> tuple[FEModel, FieldOutputController]:
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
            - The final model tree
            - The fieldoutput controller containing all processed results.
    """

    identification = "feCore"

    journal = Journal(verbose=verbose)

    job = inputfile["*job"][0]
    jobName = job["name"]

    domainSize = domainMapping[job["domain"]]

    journal.printSeperationLine()

    journal.message(
        "Setting up finite element model",
        identification,
        0,
    )

    jobInfo = dict()

    tic = getCurrentTime()
    model = FEModel(domainSize)
    model = fillFEModelFromInputFile(model, inputfile, journal)
    model.prepareYourself(journal)
    model.advanceToTime(job.get("startTime", 0.0))
    toc = getCurrentTime()
    jobInfo["model setup time"] = toc - tic

    journal.printTable(
        [
            ("Model setup time ", "{:10.4f}s".format(jobInfo["model setup time"])),
        ],
        identification,
        level=0,
    )

    printPrettyModelSummary(model, journal)
    journal.printSeperationLine()

    jobInfo["computationTime"] = 0.0

    jobInfo.update(job)
    jobInfo = loadConfiguration(jobInfo)
    for updateConfig in inputfile["*updateConfiguration"]:
        updateConfiguration(updateConfig, jobInfo, journal)

    # Create the default entries 'U' (flux) and 'P' (effort)
    for nodeField in model.nodeFields.values():
        nodeField.createFieldValueEntry("U")
        nodeField.createFieldValueEntry("P")

    model._linkFieldVariableObjects(model.nodeSets["all"])

    plotter = createPlotterFromInputFile(inputfile, journal)
    stepManager = createStepManagerFromInputFile(inputfile)
    fieldOutputController = createFieldOutputFromInputFile(inputfile, model, journal)
    fieldOutputController.initializeJob()

    outputManagers = createOutputManagersFromInputFile(
        inputfile, jobName, model, fieldOutputController, journal, plotter
    )
    for outputManager in outputManagers:
        outputManager.initializeJob()

    solvers = createSolversFromInputFile(inputfile, jobInfo, journal)

    if not solvers:
        from warnings import warn

        warn(
            "Warning, not defining a Solver is deprecated; Define solver using *solver keyword",
            DeprecationWarning,
            stacklevel=2,
        )

    defaultSolver = getSolverByName(job.get("solver", "NIST"))
    solvers["default"] = defaultSolver(jobInfo, journal)

    try:
        for step in stepManager.dequeueStep(jobInfo, model, fieldOutputController, journal, solvers, outputManagers):
            tic = getCurrentTime()
            step.solve()
            toc = getCurrentTime()
            stepTime = toc - tic
            jobInfo["computationTime"] += stepTime

            journal.printTable(
                [
                    ("Step computation time", "{:10.4f}s".format(stepTime)),
                ],
                identification,
                level=0,
            )

    except KeyboardInterrupt:
        print("")
        journal.errorMessage("Interrupted by user", identification)

    except StepFailed:
        print("")
        journal.errorMessage("Simulation failed", identification)

    except Exception as e:
        print("")
        journal.errorMessage("Simulation failed due to unhandled exception", identification)
        raise e

    finally:
        journal.printTable(
            [
                (
                    "Job computation time",
                    "{:10.4f}s".format(jobInfo["computationTime"]),
                ),
            ],
            identification,
            level=0,
            printHeaderRow=False,
        )

        fieldOutputController.finalizeJob()
        for manager in outputManagers:
            manager.finalizeJob()

        plotter.finalize()
        if not suppressPlots:
            plotter.show()

    return model, fieldOutputController
