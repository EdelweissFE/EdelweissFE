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
Created on Sat May  6 20:05:41 2017

@author: Matthias Neuner
"""

from abc import ABC, abstractmethod


class StepActionBase(ABC):
    """This is the abstract base class for all step actions managers.
    User defined step actions must implement the abstract methods."""

    identification = "StepActionBase"

    @abstractmethod
    def __init__(self, name, definition, jobInfo, modelInfo, fieldOutputController, journal):
        pass

    @abstractmethod
    def updateStepAction(self, definition):
        """is called when an updated definition is present for a new step"""
        pass

    @abstractmethod
    def finishStep(self, U, P):
        """is called when a step successfully finished"""
        pass
