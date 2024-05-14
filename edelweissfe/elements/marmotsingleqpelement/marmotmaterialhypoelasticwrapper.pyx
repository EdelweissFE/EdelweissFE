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

import numpy as np

cimport cython
cimport libcpp.cast
cimport numpy as np

from edelweissfe.elements.marmotsingleqpelement.marmot cimport (
    MarmotMaterial,
    MarmotMaterialHypoElastic,
    StateView,
    createMaterial,
)

from edelweissfe.utils.exceptions import CutbackRequest

from libc.stdlib cimport free, malloc
from libcpp.memory cimport allocator, make_unique, unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef class MarmotMaterialHypoElasticWrapper:

    cdef MarmotMaterialHypoElastic* marmotMaterialHypoElastic

    cdef readonly list fields
    cdef readonly int nU

    cdef double[::1] stateVars, stressInStateVars, strainInStateVars, stateVarsMaterial, materialProperties, dStress_dStrainInStateVars

    def __init__(self, ):

        self.fields = ["strain symmetric"]
        self.nU = 6


    def createMaterial(self, materialName, materialProperties):

        cdef MarmotMaterial* marmotMaterial = createMaterial(materialName, materialProperties)

        self.marmotMaterialHypoElastic = <MarmotMaterialHypoElastic*> (marmotMaterial)

        if self.marmotMaterialHypoElastic == NULL:
            raise ValueError("Casting to MarmotMaterialHypoElastic failed!")

    def initializeYourself(self):
        self.marmotMaterialHypoElastic.initializeYourself()

    def setCharacteristicElementLength(self, double[::1] values):
        self.marmotMaterialHypoElastic.setCharacteristicElementLength(values[0])

    def computeYourself(self,
                         double[::1] Ke,
                         double[::1] Pe,
                         const double[::1] U,
                         const double[::1] dU,
                         const double[::1] time,
                         double dTime):

        cdef double pNewDT
        pNewDT = 1e36

        cdef double [::1] stress = Pe
        cdef double [::1] dStress_dStrain = Ke
        cdef const double [::1] dStrain = dU

        stress[:] = self.stressInStateVars

        self.marmotMaterialHypoElastic.computeStress(
                &stress[0],
                &dStress_dStrain[0],
                &dStrain[0],
                &time[0],
                dTime,
                pNewDT)

        cdef int i = 0
        for i in range(6):
            self.strainInStateVars[i] += dU [i]

        self.stressInStateVars[:] = stress
        self.dStress_dStrainInStateVars[:] = dStress_dStrain

        if pNewDT < 1.0:
            raise CutbackRequest("Material requests for a cutback!", pNewDT)

    def getNumberOfRequiredStateVars(self,):

        numberOfRequiredStateVarsMaterial = self.marmotMaterialHypoElastic.getNumberOfRequiredStateVars()
        numberOfRequiredStateVarsOverhead = 6 + 6 + 36
        return numberOfRequiredStateVarsOverhead + numberOfRequiredStateVarsMaterial

    def assignStateVars(self, stateVars):

        self.stateVars = stateVars
        self.stressInStateVars = self.stateVars[0:6]
        self.strainInStateVars = self.stateVars[6:12]
        self.dStress_dStrainInStateVars = self.stateVars[12:48]
        self.stateVarsMaterial = self.stateVars[48:]

        self.marmotMaterialHypoElastic.assignStateVars(&self.stateVarsMaterial[0],
                                                       len(self.stateVarsMaterial))

    def getResultArray(self, result, getPersistentView=True):

        if result == "stress":
            return np.array(self.stressInStateVars, copy= not getPersistentView)

        if result == "strain":
            return np.array(self.strainInStateVars, copy= not getPersistentView)

        if result == "dStress_dStrain":
            return np.array(self.dStress_dStrainInStateVars, copy= not getPersistentView)

        cdef string result_ =  result.encode('UTF-8')

        cdef StateView res = self.marmotMaterialHypoElastic.getStateView(result_)

        cdef double[::1] theView = <double[:res.stateSize]> ( res.stateLocation )

        return np.array(  theView, copy= not getPersistentView)

    def __dealloc__(self):

        del self.marmotMaterialHypoElastic
