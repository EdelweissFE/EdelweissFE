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
# Created on Wed Aug 31 08:35:06 2022

# @author: matthias

from libc.stdlib cimport free, malloc
from libcpp.memory cimport allocator, make_unique, unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector

from edelweissfe.elements.marmotsingleqpelement.marmot cimport (
    MarmotMaterial,
    MarmotMaterialFactory,
)


cdef MarmotMaterial* createMaterial(materialName, materialProperties) except NULL:

    cdef double[::1] matPropsView = materialProperties

    cdef MarmotMaterial* marmotMaterial

    try:
        marmotMaterial = MarmotMaterialFactory.createMaterial(
                            MarmotMaterialFactory.getMaterialCodeFromName(
                            materialName.upper().encode('UTF-8')),
                            &matPropsView[0],
                            materialProperties.shape[0],
                            0)

    except IndexError:
        raise NotImplementedError("Marmot material {:} not found in library.".format(materialName))

    return marmotMaterial
