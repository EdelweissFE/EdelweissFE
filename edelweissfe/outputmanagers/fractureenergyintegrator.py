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
# Created on Tue Dec 17 08:26:01 2019

# @author: Matthias Neuner

import numpy as np

from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from edelweissfe.utils.math import createMathExpression

"""
A simple integrator to compute the fracture energy by integrating a load-displacement curve.

.. code-block:: edelweiss
    :caption: Example:

    *output, type=fractureenergyintegrator, jobName=myJob, name=gfi
        forceFieldOutput=RF, displacementFieldOutput=U, fractureArea='100.0*1.0'
"""

documentation = {
    "forceFieldOutput": "fieldOutput for force (with time history)",
    "displacementFieldOutput": "fieldOutput for displacement (with time history)",
    "fractureArea": "(math expression for) area of fracture",
}


class OutputManager(OutputManagerBase):
    """Simple Integrator for fracture energy"""

    identification = "FEI"
    printTemplate = "{:}, {:}: {:}"

    def __init__(self, name, model, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []
        self.fieldOutputController = fieldOutputController

    def updateDefinition(self, **kwargs: dict):
        self.fpF = self.fieldOutputController.fieldOutputs[kwargs["forceFieldOutput"]]
        self.fpU = self.fieldOutputController.fieldOutputs[kwargs["displacementFieldOutput"]]
        self.A = createMathExpression(kwargs["fractureArea"])(0.0)
        self.fractureEnergy = 0.0

    def initializeJob(self):
        pass

    def initializeStep(self, step):
        pass

    def finalizeIncrement(self, **kwargs):
        pass

    def finalizeFailedIncrement(self, **kwargs):
        pass

    def finalizeStep(
        self,
    ):
        pass

    def finalizeJob(
        self,
    ):
        self.fractureEnergy = np.trapz(self.fpF.getResultHistory(), x=self.fpU.getResultHistory()) / self.A
        self.journal.message(
            "integrated fracture energy: {:3.4f}".format(self.fractureEnergy),
            self.identification,
        )
