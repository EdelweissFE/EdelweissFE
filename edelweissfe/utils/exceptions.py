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
# Created on Fri May 12 11:02:12 2017

# @author: Matthias Neuner
"""
This module contains a collection of commonly used exceptions.
"""


class WrongDomain(Exception):
    """Is thrown when a required module does not fit the spatial domain."""


class StepFailed(Exception):
    """Is thrown when the computation of a step fails."""


class CutbackRequest(Exception):
    """The evaluation of the residuel fails due to a too large time increment."""

    def __init__(self, message, cutbackSize):
        super().__init__(message)
        self.cutbackSize = float(cutbackSize)


class ReachedMaxIterations(Exception):
    """The maximum number of nonlinear iterations as attained."""


class ReachedMaxIncrements(Exception):
    """The maximum number of incremeents within a step as attained."""


class ReachedMinIncrementSize(Exception):
    """The minimum size of a incremeent within a step as attained."""


class DivergingSolution(Exception):
    """The solutions seems to be diverging within the nonlinear solving scheme."""


class ConditionalStop(Exception):
    """Simulation stops sucessfully."""


class InputExecption(Exception):
    """Invalid input file parameters."""
