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

from fe.config.outputmanagers import getOutputManagerClass
from fe.models.femodel import FEModel
from fe.config.outputmanagers import getOutputManagerClass
from fe.utils.fieldoutput import FieldOutputController
from fe.utils.plotter import Plotter
from fe.journal.journal import Journal


def createOutputManagersFromInput(
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
