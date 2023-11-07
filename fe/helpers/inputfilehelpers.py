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


def createFieldOutputFromInputFile(inputfile: dict, model: FEModel, journal: Journal) -> FieldOutputController:
    """Convenience helper function
    to create the FieldOutputController instance using the *fieldOutput keyword.

    Parameters
    ----------
    inputfile
        The processed inputfile in dictionary form.
    model
        The FEModel tree instance.
    journal
        The Journal for logging purposes.

    Returns
    -------
    FieldOutputController
        The configured FieldOutputController instance.
    """
    fieldOutputController = FieldOutputController(model, journal)
    if inputfile["*fieldOutput"]:
        for definition in inputfile["*fieldOutput"]:
            for defLine in definition["data"]:
                kwargs = convertLineToStringDictionary(defLine)
                if "elSet" in kwargs:
                    kwargs["elSet"] = model.elementSets[kwargs["elSet"]]
                if "nSet" in kwargs:
                    kwargs["nSet"] = model.nodeSets[kwargs["nSet"]]
                name = kwargs.pop("name")

                theType = kwargs.pop("create")

                if theType == "perNode":
                    field = kwargs.pop("field")
                    nodeField = model.nodeFields[field]
                    result = kwargs.pop("result")

                    subset = None
                    if "nSet" in kwargs:
                        subset = kwargs.pop("nSet")
                    elif "elSet" in kwargs:
                        subset = kwargs.pop("elSet")

                    if subset:
                        nodeField = nodeField.subset(subset)

                    fieldOutputController.addPerNodeFieldOutput(name, nodeField, result, **kwargs)

                elif theType == "perElement":
                    elSet = kwargs.pop("elSet")
                    result = kwargs.pop("result")
                    fieldOutputController.addPerElementFieldOutput(name, elSet, result, **kwargs)

                elif theType == "fromExpression":
                    # elSet = kwargs.pop("elSet")
                    # result = kwargs.pop("result")
                    fieldOutputController.addExpressionFieldOutput(name, **kwargs)

                else:
                    raise Exception("Invalid FieldOuput request: {:}".format(theType))

    return fieldOutputController


def fillFEModelFromInputFile(model: FEModel, inputfile: dict, journal: Journal) -> FEModel:
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


def createStepManagerFromInputFile(inputfile: dict):
    """Convenience helper function
    to create and fill a StepManager using *step definitions from the input file.

    Parameters
    ----------
    inputfile
        The processed inputfile in dictionary form.

    Returns
    -------
    StepManager
        The StepManager with enqueued StepDefinitions.
    """
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


def createSolversFromInputFile(inputfile: dict, jobInfo: dict, journal: Journal) -> dict:
    """Convenience helper function
    to create instances of nonlinear solvers using the *solver keyword.

    Parameters
    ----------
    inputfile
        The processed inputfile in dictionary form.
    jobInfo
        Additional informations about the job
    journal
        The Journal instance for logging purposes.

    Returns
    -------
    dict
        The dictionary containing the solver instances.
    """
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
    """Convenience helper function
    to create output managers using the *output keyword.

    Parameters
    ----------
    inputfile
        The processed inputfile in dictionary form.
    defaultName
        The default name prefix for output managers without specified name.
    fieldOutputController
        The FieldOutputController instance.
    journal
        The Journal instance for logging purposes.
    plotter
        The plotter instance.

    Returns
    -------
    list
        The list containing the OutputManager instances.
    """
    jobName = defaultName
    outputManagers = []

    for outputDef in inputfile["*output"]:
        OutputManager = getOutputManagerClass(outputDef["type"].lower())
        managerName = outputDef.get("name", defaultName + outputDef["type"])
        definitionLines = outputDef["data"]

        outputManager = OutputManager(managerName, model, fieldOutputController, journal, plotter)

        for defLine in definitionLines:
            kwargs = convertLineToStringDictionary(defLine)
            if "elSet" in kwargs:
                kwargs["elSet"] = model.elementSets[kwargs["elSet"]]
            if "nSet" in kwargs:
                kwargs["nSet"] = model.nodeSets[kwargs["nSet"]]
            if "fieldOutput" in kwargs:
                kwargs["fieldOutput"] = fieldOutputController.fieldOutputs[kwargs["fieldOutput"]]

            outputManager.updateDefinition(**kwargs)

        outputManagers.append(outputManager)

    return outputManagers


def createPlotterFromInputFile(inputfile: dict, journal: Journal) -> Plotter:
    """Convenience helper function
    to create and configure the Plotter instance using the *configurePlots and *exportPlots keywords

    Parameters
    ----------
    inputfile
        The processed inputfile in dictionary form.
    journal
        The Journal instance for logging purposes.

    Returns
    -------
    Plotter
        The resulting plotter instance
    """
    plotConfigurations = [
        convertLineToStringDictionary(c) for configEntry in inputfile["*configurePlots"] for c in configEntry["data"]
    ]

    exportJobs = [
        convertLineToStringDictionary(c) for configEntry in inputfile["*exportPlots"] for c in configEntry["data"]
    ]

    plotter = Plotter(journal, plotConfigurations, exportJobs)

    return plotter
