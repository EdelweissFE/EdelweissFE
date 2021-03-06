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
Created on Fri Feb 17 17:09:35 2017

@author: Matthias Neuner
"""

from abc import ABC, abstractmethod


class OutputManagerBase(ABC):
    """This is the abstract base class for all output managers.
    User defined output managers must implement the abstract methods."""

    identification = "OutputManagerBase"

    @abstractmethod
    def __init__(self, name, definitionLines, jobInfo, modelInfo, fieldOutputController, journal, plotter):
        pass

    @abstractmethod
    def initializeStep(self, step, stepActions, stepOptions):
        pass

    @abstractmethod
    def finalizeIncrement(self, U, P, increment):
        pass

    @abstractmethod
    def finalizeStep(self, U, P, time):
        pass

    @abstractmethod
    def finalizeJob(self, U, P):
        pass
