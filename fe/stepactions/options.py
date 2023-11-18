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
# Created on Mon Jan 23 13:03:09 2017

# @author: Matthias Neuner
"""This stepaction serves as a simple case insensitive container for
storing step options for various modules.
"""

from fe.stepactions.base.stepactionbase import StepActionBase
import numpy as np
import sympy as sp
from collections import defaultdict


class StepAction(StepActionBase):
    def __init__(self, name, options, jobInfo, model, fieldOutputController, journal):
        self.name = name
        self.options = options
        self.updateStepAction(options, jobInfo, model, fieldOutputController, journal)

    def __contains__(self, key):
        """We may work with this action like a dictionary."""

        return key in self.options

    def __getitem__(self, key):
        """We may work with this action like a dictionary."""
        return self.options[key]

    def __setitem__(self, key, val):
        """We may work with this action like a dictionary."""

        self.options[key] = val

    def get(self, key, default):
        """We may work with this action like a  dictionary."""

        return self.options.get(key, default)

    def updateStepAction(self, options, jobInfo, model, fieldOutputController, journal):
        self.options.update(options)
