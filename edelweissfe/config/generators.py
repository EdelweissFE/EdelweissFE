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
# Created on Wed Apr 12 15:37:28 2017

# @author: Matthias Neuner
"""
Mesh generators are used to create meshes,
using different methods.

Keyword: ``*generator``

.. code-block:: edelweiss
    :caption: Example:

    *modelGenerator, generator=theGeneratorType, name=myGeneratorName
        multiple lines of defintion ...
        multiple lines of defintion ...
        multiple lines of defintion ...
"""

import importlib


def getGeneratorFunction(name: str) -> type:
    """Get the function type of the requested generator.

    Parameters
    ----------
    name
        The name of the generator class type to load.
    Returns
    -------
    type
        The generator function type.
    """

    module = importlib.import_module("edelweissfe.generators." + name.lower())
    return module.generateModelData
