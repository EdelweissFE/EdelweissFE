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
#  Paul Hofer Paul.Hofer@uibk.ac.at
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

"""This stepaction serves as a simple case insensitive container for
storing step options for various modules.
"""

from edelweissfe.stepactions.base.stepactionbase import StepActionBase
from edelweissfe.utils.caseinsensitivedict import CaseInsensitiveDict


class StepAction(StepActionBase):
    def __init__(self, name, options, jobInfo, model, fieldOutputController, journal):
        self.name = name
        self.options = CaseInsensitiveDict(options)
        self.updateStepAction(options, jobInfo, model, fieldOutputController, journal)

    def updateStepAction(self, options, jobInfo, model, fieldOutputController, journal):
        self.options.update(options)

    def __contains__(self, *args):
        """wrapper method for CaseInsensitiveDict"""
        return self.options.__contains__(*args)

    def __getitem__(self, *args):
        """wrapper method for CaseInsensitiveDict"""
        self.options.__getitem__(*args)

    def __setitem__(self, *args):
        """wrapper method for CaseInsensitiveDict"""
        self.options.__setitem__(*args)

    def get(self, *args):
        """wrapper method for CaseInsensitiveDict"""
        return self.options.get(*args)
