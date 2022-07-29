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
# Created on Sat Jul 22 21:26:01 2017

# @author: Matthias Neuner
"""
A conditional stop conditions wenn an expression becomes true.

Datalines:
"""
documentation = {
    "fieldOutput": "fieldOutput to be monitored",
    "f(x)": "(optional), apply a model accesible expression on the fieldOutput",
}

from fe.outputmanagers.outputmanagerbase import OutputManagerBase
from fe.utils.misc import stringDict
from fe.utils.math import createMathExpression
from fe.utils.exceptions import ConditionalStop


class OutputManager(OutputManagerBase):
    identification = "ConditionalStop"
    printTemplate = "{:}, {:}: {:}"

    def __init__(self, name, definitionLines, jobInfo, modelInfo, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []
        self.fieldOutputController = fieldOutputController

        for defline in definitionLines:
            entry = {}
            defDict = stringDict(defline)
            entry["fieldOutput"] = fieldOutputController.fieldOutputs[defDict["fieldOutput"]]
            entry["f(x)"] = createMathExpression(defDict.get("f(x)", "x"))
            self.monitorJobs.append(entry)

    def initializeStep(self, step, stepActions, stepOptions):
        pass

    def finalizeIncrement(self, U, P, increment):
        for nJob in self.monitorJobs:
            result = nJob["f(x)"](nJob["fieldOutput"].getResultHistory())

            if result == True:
                raise ConditionalStop()

    def finalizeStep(self, U, P, time):
        pass

    def finalizeJob(self, U, P):
        pass
