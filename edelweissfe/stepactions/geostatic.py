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
# Created on Tue May  9 19:52:53 2017

# @author: Matthias Neuner

# @ Magdalena
"""
Initialize materials to an geostatic stress state
"""

import numpy as np

from edelweissfe.stepactions.base.stepactionbase import StepActionBase

documentation = {
    "p1": "sig_x=sig_y=sig_z in first point",
    "h1": "y coordinate of first point, default=1.0",
    "p2": "s11=s22=s33 in second point, default=p1",
    "h2": "y coordinate of second point, default=-1.0",
    "xLateral": "ratio of sig_x/sig_y, default=1.0",
    "zLateral": "ratio of sig_z/sig_y, default=1.0",
}


class StepAction(StepActionBase):
    """Initializes elements of set with an Abaqus-like geostatic stress state.
    Is automatically deactivated at the end of the step."""

    def __init__(self, name, action, jobInfo, model, fieldOutputController, journal):
        self.name = name

        self.geostaticElements = model.elementSets[action.get("elSet", "all")]
        self.p1 = float(action["p1"])
        self.p2 = float(action.get("p2", self.p1))
        self.level1 = float(action.get("h1", 1.0))
        self.level2 = float(action.get("h2", -1.0))
        self.xLateral = float(action.get("xLateral", 1.0))
        self.zLateral = float(action.get("zLateral", 1.0))

        self.geostaticDefinition = np.array(
            [
                self.p1,
                self.level1,
                self.p2,
                self.level2,
                self.xLateral,
                self.zLateral,
            ]
        )

        self.journal = journal

        self.active = True

    def applyAtStepEnd(self, model, stepMagnitude=None):
        if not self.active:
            return

        self.journal.printSeperationLine()
        self.journal.message("End of geostatic step -- displacements are reset", self.name)
        self.journal.printSeperationLine()

        model.nodeFields["displacement"]["U"][:] = 0
        # U[self.theDofManager.indicesOfFieldsInDofVector["displacement"]] = 0.0

        self.active = False

    def applyAtIterationStart(
        self,
    ):
        if not self.active:
            return

        for el in self.geostaticElements:
            el.setInitialCondition("geostatic stress", self.geostaticDefinition)
