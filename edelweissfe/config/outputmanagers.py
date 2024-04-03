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
# Created on Thu Jan 26 09:27:14 2017

# @author: Matthias Neuner
"""
Outputmanagers are used for monitoring and exporting results to variuous formats.
They are defined globally:

.. code-block:: edelweiss

    *output, type=ensight, jobName=cps4job, name=myEnsightExport
        create=perNode,    elSet=all, fieldOutput=Displacement
        create=perElement, elSet=all, fieldOutput=Stress
        create=perElement, elSet=all, fieldOutput=Strain

"""
import importlib


def getOutputManagerClass(name: str) -> type:
    """Get the class type of the requested output manager.

    Parameters
    ----------
    name
        The name of the output manager to load.

    Returns
    -------
    type
        The output manager class type.
    """

    module = importlib.import_module("edelweissfe.outputmanagers." + name.lower())
    return module.OutputManager
