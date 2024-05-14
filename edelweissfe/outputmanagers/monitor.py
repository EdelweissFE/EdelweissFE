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

from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from edelweissfe.utils.math import createMathExpression

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


class OutputManager(OutputManagerBase):
    """Simple monitor for nodes, nodeSets, elements and elementSets"""

    identification = "Monitor"
    printTemplate = "{:}: {:}"

    def __init__(self, name, model, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []
        self.fieldOutputController = fieldOutputController

    def updateDefinition(self, **kwargs: dict):
        entry = dict()
        entry["fieldOutput"] = kwargs["fieldOutput"]
        entry["f(x)"] = createMathExpression(kwargs.get("f(x)", "x"))
        self.monitorJobs.append(entry)

    def initializeJob(self):
        pass

    def initializeStep(self, step):
        pass

    def finalizeIncrement(self, **kwargs):
        for nJob in self.monitorJobs:
            result = nJob["f(x)"](nJob["fieldOutput"].getLastResult())
            self.journal.message(
                self.printTemplate.format(nJob["fieldOutput"].name, result),
                self.identification,
            )

    def finalizeFailedIncrement(self, **kwargs):
        pass

    def finalizeStep(
        self,
    ):
        pass

    def finalizeJob(
        self,
    ):
        pass
