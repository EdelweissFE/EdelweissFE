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
Created on Wed Sep 13 08:50:46 2017

@author: Matthias Neuner
"""

import numpy as np


def sigPrinc(x):
    return np.linalg.eig([[x[0], x[3], x[4]], [x[3], x[1], x[5]], [x[4], x[5], x[2]]])[0]


mathModules = {
    "sigPrinc": sigPrinc,
    "mean": np.mean,
    "max": np.max,
    "min": np.min,
    "abs": np.abs,
    "sum": np.sum,
    "linalg": np.linalg,
}


def createFunction(expression, *argnames, **kwargs):
    """Create a function from a string expression"""
    scope = {**locals(), **kwargs}
    return lambda *args: eval(expression, globals(), {**dict(zip(argnames, args)), **kwargs})


def createModelAccessibleFunction(expression, model, *argnames, **kwargs):
    """Create a function from a string expression, which can access the complete model any given objects"""
    kwargs = {**kwargs, "model": model}
    return createFunction(expression, *argnames, **kwargs)


def evalExpression(expression, **kwargs):
    """Create a function from a string expression"""
    return eval(expression, globals(), kwargs)


def evalModelAccessibleExpression(expression, model, *args, **kwargs):
    """Evalualate a string expression, which can access the complete model"""
    kwargs = {**kwargs, "model": model}
    return evalExpression(expression, **kwargs)


def execExpression(expression, **kwargs):
    """Create a function from a string expression"""
    return exec(expression, globals(), kwargs)


def execModelAccessibleExpression(expression, model, *args, **kwargs):
    """Evalualate a string expression, which can access the complete model"""
    kwargs = {**kwargs, "model": model}
    return execExpression(expression, **kwargs)


def createMathExpression(expression, symbol="x"):
    return lambda x: eval(expression, globals(), {symbol: x, **mathModules})
