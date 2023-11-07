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
Useful, e.g., for indirect displacement control.

.. code-block:: edelweiss
    :caption: Example: 

    *output, type=conditionalstop, jobName=myJob, name=condStop
        stop='fieldOutputs["damage"]  >= .99'
        stop='fieldOutputs["displacement"]  < -5'
"""
documentation = {"stop": "model accessible function describing the stop condition"}

from fe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from fe.utils.misc import convertLineToStringDictionary
from fe.utils.math import createModelAccessibleFunction
from fe.utils.exceptions import ConditionalStop


class OutputManager(OutputManagerBase):
    identification = "ConditionalStop"
    printTemplate = "{:}, {:}: {:}"

    def __init__(self, name, model, fieldOutputController, journal, plotter):
        self.model = model
        self.journal = journal
        self.monitorJobs = []
        self.fieldOutputController = fieldOutputController

    def updateDefinition(self, **kwargs: dict):
        # for defline in definitionLines:
        entry = {}
        entry["stop"] = createModelAccessibleFunction(
            kwargs["stop"], model, fieldOutputs=fieldOutputController.fieldOutputs
        )
        self.monitorJobs.append(entry)

    def initializeJob(self):
        pass

    def initializeStep(self, step):
        pass

    def finalizeIncrement(self, **kwargs):
        for nJob in self.monitorJobs:
            if nJob["stop"]():
                raise ConditionalStop()

    def finalizeFailedIncrement(self, **kwargs):
        pass

    def finalizeStep(self):
        pass

    def finalizeJob(self):
        pass
