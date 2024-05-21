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

from edelweissfe.config.generators import getGeneratorFunction
from edelweissfe.config.outputmanagers import getOutputManagerClass
from edelweissfe.config.solvers import getSolverByName
from edelweissfe.generators.abqmodelconstructor import AbqModelConstructor
from edelweissfe.journal.journal import Journal
from edelweissfe.models.femodel import FEModel
from edelweissfe.steps.stepmanager import (
    StepActionDefinition,
    StepDefinition,
    StepManager,
)
from edelweissfe.utils.fieldoutput import FieldOutputController
from edelweissfe.utils.math import createMathExpression, createModelAccessibleFunction
from edelweissfe.utils.misc import (
    convertAssignmentsToStringDictionary,
    convertLinesToStringDictionary,
    convertLineToStringDictionary,
    isInteger,
    splitLineAtCommas,
    strToRange,
)
from edelweissfe.utils.plotter import Plotter


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

                f_of_x = kwargs.pop("f(x)", None)
                if f_of_x:
                    f_of_x = createMathExpression(f_of_x)

                f_export_of_x = kwargs.pop("f_export(x)", None)
                if f_export_of_x:
                    f_export_of_x = createMathExpression(f_export_of_x)

                saveHistory = kwargs.pop("saveHistory", True)
                if saveHistory:
                    saveHistory = bool(saveHistory)

                export = kwargs.pop("export", False)

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

                    fieldOutputController.addPerNodeFieldOutput(
                        name,
                        nodeField,
                        result,
                        saveHistory=saveHistory,
                        f_x=f_of_x,
                        export=export,
                        fExport_x=f_export_of_x,
                        **kwargs,
                    )

                elif theType == "perElement":
                    elSet = kwargs.pop("elSet")
                    result = kwargs.pop("result")

                    qp = kwargs.pop("quadraturePoint")
                    quadraturePoints = strToRange(qp) if not isInteger(qp) else [int(qp)]

                    fieldOutputController.addPerElementFieldOutput(
                        name,
                        elSet,
                        result,
                        saveHistory=saveHistory,
                        f_x=f_of_x,
                        export=export,
                        fExport_x=f_export_of_x,
                        quadraturePoints=quadraturePoints,
                    )

                elif theType == "fromExpression":

                    if "nSet" in kwargs:
                        associatedSet = kwargs.pop("nSet")
                    elif "elSet" in kwargs:
                        associatedSet = kwargs.pop("elSet")
                    else:
                        raise Exception("All FieldOuputs must be associated with a set!")

                    theExpression = createModelAccessibleFunction(kwargs["expression"], model)

                    fieldOutputController.addExpressionFieldOutput(
                        associatedSet,
                        theExpression,
                        name,
                        saveHistory,
                        f_of_x,
                        export=export,
                        fExport_x=f_export_of_x,
                    )
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
        try:
            solverName = solverDefinition["name"]
        except KeyError:
            raise KeyError("Solver definition missing name")
        try:
            solverType = solverDefinition["solver"]
        except KeyError:
            raise KeyError(f"Missing type definition for solver {solverName}. Specify solver type with solver=...")

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
