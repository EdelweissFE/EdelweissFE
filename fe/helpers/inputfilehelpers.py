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

from fe.utils.fieldoutput import FieldOutputController
from fe.utils.misc import convertLineToStringDictionary
from fe.models.femodel import FEModel
from fe.journal.journal import Journal
from fe.config.generators import getGeneratorFunction
from fe.utils.abqmodelconstructor import AbqModelConstructor
from fe.variables.scalarvariable import ScalarVariable
from fe.utils.misc import convertAssignmentsToStringDictionary, splitLineAtCommas, convertLinesToStringDictionary
from fe.steps.stepmanager import StepManager, StepActionDefinition, StepDefinition
from fe.utils.misc import convertLinesToStringDictionary
from fe.config.solvers import getSolverByName
from fe.config.outputmanagers import getOutputManagerClass
from fe.config.outputmanagers import getOutputManagerClass
from fe.utils.fieldoutput import FieldOutputController
from fe.utils.plotter import Plotter


def createFieldOutputFromInputFile(inputfile: dict, model, journal):
    fieldOutputController = FieldOutputController()
    if inputfile["*fieldOutput"]:
        definition = inputfile["*fieldOutput"][0]
        for defLine in definition["data"]:
            fpDef = convertLineToStringDictionary(defLine)
            name = fpDef.pop("name")
            fieldOutputController.addFieldOutput(name, model, journal, **fpDef)

    return fieldOutputController


def fillFEModelFromInputFile(model: FEModel, inputfile: dict, journal: Journal):
    """Convenience helper function
    to fill an existing (possibly empty) FEModel using the input file and generators.

    Parameters
    ----------
    FEModel
        The model tree to be filled.
    input
        The processed inputfile in dictionary form.
    journal
        The Journal for logging purposes.

    Returns
    -------
    FEModel
        The updated, filled model tree.
    """

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

    return model


def createStepManagerFromInputFile(inputfile):
    stepManager = StepManager()

    for stepLine in inputfile["*step"]:
        stepType = stepLine.pop("type", "AdaptiveStep")
        stepActionLines = stepLine.pop("data")

        stepActionDefinitions = []

        for line in stepActionLines:
            module, *definition = splitLineAtCommas(line)
            kwargs = convertAssignmentsToStringDictionary(definition)

            name = kwargs.pop("name", kwargs.pop("category", module))

            stepActionDefinitions.append(StepActionDefinition(name, module, kwargs))

        stepDefinition = StepDefinition(stepType, stepLine, stepActionDefinitions)

        stepManager.enqueueStepDefinition(stepDefinition)

    return stepManager


def createSolversFromInputFile(inputfile: dict, jobInfo: dict, journal: Journal):
    solvers = {}
    for solverDefinition in inputfile["*solver"]:
        solverType = solverDefinition["solver"]
        solverName = solverDefinition["name"]
        solverData = solverDefinition["data"]

        Solver = getSolverByName(solverType)

        solverData = convertLinesToStringDictionary(solverData)

        solvers[solverName] = Solver(jobInfo, journal, **solverData)

    return solvers


def createOutputManagersFromInputFile(
    inputfile: dict,
    defaultName: str,
    model: FEModel,
    fieldOutputController: FieldOutputController,
    journal: Journal,
    plotter: Plotter,
) -> list:
    jobName = defaultName
    outputmanagers = []

    for outputDef in inputfile["*output"]:
        OutputManager = getOutputManagerClass(outputDef["type"].lower())
        managerName = outputDef.get("name", defaultName + outputDef["type"])
        definitionLines = outputDef["data"]
        outputmanagers.append(
            OutputManager(managerName, definitionLines, model, fieldOutputController, journal, plotter)
        )

    return outputmanagers


def createPlotterFromInputFile(inputfile: dict, journal: Journal):
    plotConfigurations = [
        convertLineToStringDictionary(c) for configEntry in inputfile["*configurePlots"] for c in configEntry["data"]
    ]

    exportJobs = [
        convertLineToStringDictionary(c) for configEntry in inputfile["*exportPlots"] for c in configEntry["data"]
    ]

    plotter = Plotter(journal, plotConfigurations, exportJobs)

    return plotter
