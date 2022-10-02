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
"""
Created on Thu Apr 13 14:08:32 2017

@author: Matthias Neuner
"""


from fe.outputmanagers.base.outputmanagerbase import OutputManagerBase

from fe.utils.misc import convertLinesToStringDictionary
import numpy as np


class OutputManager(OutputManagerBase):
    identification = "TimeMonitor"

    def __init__(self, name, definitionLines, jobInfo, model, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []
        defDict = convertLinesToStringDictionary(definitionLines)
        self.exportFile = defDict["export"]
        self.timeVals = []

    def initializeSimulation(self, model):
        pass

    def initializeStep(self, step, stepActions):
        pass

    def finalizeIncrement(self, U, P, increment, **kwargs):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        self.timeVals.append(totalTime + dT)

    def finalizeFailedIncrement(self, **kwargs):
        pass

    def finalizeStep(self, U, P, time):
        pass

    def finalizeJob(
        self,
        U,
        P,
    ):
        np.savetxt("{:}.csv".format(self.exportFile), np.asarray(self.timeVals).T)
