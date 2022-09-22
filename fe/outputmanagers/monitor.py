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
A simple monitor to observe results (fieldOutputs) in the console during analysis.

.. code-block:: edelweiss
    :caption: Example:

    *output, type=monitor, jobName=cpe4job, name=omegaMon
        fieldOutput=omega, f(x)='max(x)'
"""

documentation = {
    "fieldOutput": "Name of the field output to be monitored",
    "f(x)": "(Optional), apply a model accessible function on the result",
}

from fe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from fe.utils.misc import convertLineToStringDictionary
from fe.utils.math import createMathExpression


class OutputManager(OutputManagerBase):
    """Simple monitor for nodes, nodeSets, elements and elementSets"""

    identification = "Monitor"
    printTemplate = "{:}, {:}: {:}"

    def __init__(self, name, definitionLines, jobInfo, model, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []
        self.fieldOutputController = fieldOutputController

        for defline in definitionLines:
            entry = {}
            defDict = convertLineToStringDictionary(defline)
            entry["fieldOutput"] = fieldOutputController.fieldOutputs[defDict["fieldOutput"]]
            entry["f(x)"] = createMathExpression(defDict.get("f(x)", "x"))
            self.monitorJobs.append(entry)

    def initializeSimulation(self, model):
        pass

    def initializeStep(self, step, stepActions):
        pass

    def finalizeIncrement(self, U, P, increment, **kwargs):
        for nJob in self.monitorJobs:
            result = nJob["f(x)"](nJob["fieldOutput"].getLastResult())
            self.journal.message(
                self.printTemplate.format(nJob["fieldOutput"].name, nJob["fieldOutput"].type, result),
                self.identification,
            )

    def finalizeFailedIncrement(self, **kwargs):
        pass

    def finalizeStep(self, U, P, time):
        pass

    def finalizeJob(self, U, P):
        pass
