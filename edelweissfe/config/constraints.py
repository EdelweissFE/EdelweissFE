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
# Created on Sun May 26 11:27:14 2017

# @author: Matthias Neuner
"""
Constraints are defined globally,
and may be used for introducing additional scalar
equations into the global equation system.

Example syntax:

.. code-block:: edelweiss

    *constraint, type=rigidbody, name=myRigigBody
        nSet=right
        referencePoint=rBottom
"""
import importlib


def getConstraintClass(name: str) -> type:
    """Get the class type of the requested constraint.

    Parameters
    ----------
    name
        The name of the constraint to load.

    Returns
    -------
    type
        The constraint class type.
    """

    module = importlib.import_module("edelweissfe.constraints." + name.lower())
    return module.Constraint
