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
Created on Mon May  1 19:38:14 2017

@author: Matthias Neuner
"""
import importlib
import textwrap


def printDocumentation(module):
    try:
        mod = importlib.import_module("edelweissfe." + module)
        print(mod.__doc__)
        wrapper = textwrap.TextWrapper(width=80, replace_whitespace=False)
        for key, val in mod.documentation.items():
            wrapper.initial_indent = "    {:26}".format(key)
            wrapper.subsequent_indent = " " * len(wrapper.initial_indent)
            print(wrapper.fill(val))

    except ModuleNotFoundError as e:
        print(e)
    except AttributeError as e:
        print("documentation for module {:} not found:".format(module), str(e))
