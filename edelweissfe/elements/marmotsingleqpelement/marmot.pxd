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

from libcpp.string cimport string
from libcpp.vector cimport vector


cdef extern from "Marmot/Marmot.h" namespace "MarmotLibrary" nogil:

    cdef cppclass MarmotMaterialFactory:
        @staticmethod
        int getMaterialCodeFromName(const string& materialName) except +IndexError

        @staticmethod
        MarmotMaterial* createMaterial(int materialCode, const double* materialProperties, int nMaterialProperties, int materialNumber) except +IndexError

cdef extern from "Marmot/MarmotUtils.h":
    cdef struct StateView:
        double *stateLocation
        int stateSize

cdef extern from "Marmot/MarmotMaterial.h":
    cdef cppclass MarmotMaterial nogil:
        pass

cdef extern from "Marmot/MarmotMaterialHypoElastic.h":
    cdef cppclass MarmotMaterialHypoElastic nogil:

        void assignStateVars( double* stateVars, int nStateVars )

        StateView getStateView( const string& stateName ) except +ValueError

        void initializeYourself()

        void setCharacteristicElementLength(double length)

        int getNumberOfAssignedStateVars()

        int getNumberOfRequiredStateVars()

        double* getAssignedStateVars()

        void computeStress( double*       stress,
                              double*       dStressDDStrain,
                              const double* dStrain,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) except +ValueError

cdef extern from "Marmot/MarmotMaterialGradientEnhancedHypoElastic.h":
    cdef cppclass MarmotMaterialGradientEnhancedHypoElastic nogil:

        void assignStateVars( double* stateVars, int nStateVars )

        StateView getStateView( const string& stateName ) except +ValueError

        void initializeYourself()

        void setCharacteristicElementLength(double length)

        int getNumberOfAssignedStateVars()

        int getNumberOfRequiredStateVars()

        double* getAssignedStateVars()

        void computeStress(   double*       stress,
                              double&       K_local,
                              double&       nonLocalRadius,
                              double*       dStressDDStrain,
                              double*       dK_localDDStrain,
                              double*       dStressDK,
                              const double* dStrain,
                              double        KOld,
                              double        dK,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT
                           ) except +ValueError

cdef MarmotMaterial* createMaterial(materialName, materialProperties) except NULL
