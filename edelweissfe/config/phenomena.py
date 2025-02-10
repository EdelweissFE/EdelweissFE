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
# Created on Tue Jan  10 19:10:42 2017

# @author: Matthias Neuner
"""
This module contains the central definition of
the known physics,
i.e.,

 * the available physical fields
 * their tensorial order
 * their default residual tolerances

If a new physical field should be added to EdelweissFE,
then this should happen (only) here.
"""

# field                  field type
phenomena = {
    "displacement": "vector",
    "rotation": "rotation vector",
    "micro rotation": "rotation vector",
    "thermal": "scalar",
    "nonlocal damage": "scalar",
    "nonlocal damage 2": "scalar",
    "concentration": "scalar",
    "chemical potential": "scalar",
    "strain symmetric": "symmetric tensor second order",
}


# field                  tolerance
fieldCorrectionTolerance = {
    "displacement": 1e-8,
    "rotation": 1e-3,
    "micro rotation": 1e-8,
    "nonlocal damage": 1e-8,
    "nonlocal damage 2": 1e-8,
    "concentration": 1e-1,
    "chemical potential": 1e-1,
    "strain symmetric": 1e-8,
    "scalar variables": 1e-3,
}

fluxResidualTolerance = {
    "displacement": 1e-8,
    "rotation": 1e-4,
    "micro rotation": 1e-8,
    "nonlocal damage": 1e-8,
    "nonlocal damage 2": 1e-8,
    "concentration": 1e-1,
    "chemical potential": 1e-2,
    "strain symmetric": 1e-8,
    "scalar variables": 1e-8,
}

fluxResidualToleranceAlternative = {
    "displacement": 5e-3,
    "rotation": 5e-3,
    "micro rotation": 5e-3,
    "nonlocal damage": 5e-3,
    "nonlocal damage 2": 5e-3,
    "concentration": 5e-2,
    "chemical potential": 5e-2,
    "strain symmetric": 5e-3,
    "scalar variables": 1e-8,
}

# domain                 dimensions
domainMapping = {
    "1d": 1,
    "2d": 2,
    "3d": 3,
    "axisymmetric": 2,
}


def getFieldSize(field, domainSize):
    fType = phenomena[field]
    if fType == "scalar":
        return 1
    if fType == "vector":
        return domainSize
    if fType == "rotation vector":
        if domainSize == 2:
            return 1
        elif domainSize == 3:
            return 3
    if fType == "symmetric tensor second order":
        if domainSize == 2:
            return 3
        elif domainSize == 3:
            return 6

    raise NotImplementedError("Invalid physical field {:} requested".format(field))
