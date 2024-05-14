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
    MarmotMaterialGradientEnhancedHypoElastic,
    StateView,
    createMaterial,
)

from edelweissfe.utils.exceptions import CutbackRequest

from libc.stdlib cimport free, malloc
from libcpp.memory cimport allocator, make_unique, unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef class MarmotMaterialGradientEnhancedHypoElasticWrapper:

    cdef MarmotMaterialGradientEnhancedHypoElastic* marmotMaterialGradientEnhancedHypoElastic

    cdef readonly list fields

    cdef readonly int nU

    cdef double[::1] stateVars, stressLikeInStateVars, strainLikeInStateVars, stateVarsMaterial, materialProperties, algorithmicTangentInStateVars

    def __init__(self, ):

        self.fields = ["strain symmetric", "nonlocal damage"]
        self.nU = 7

    def createMaterial(self, materialName, materialProperties):

        cdef MarmotMaterial* marmotMaterial = createMaterial(materialName, materialProperties)

        self.marmotMaterialGradientEnhancedHypoElastic = <MarmotMaterialGradientEnhancedHypoElastic*> (marmotMaterial)

        if self.marmotMaterialGradientEnhancedHypoElastic == NULL:
            raise ValueError("Casting to MarmotMaterialGradientEnhancedHypoElastic failed!")

    def initializeYourself(self):
        self.marmotMaterialGradientEnhancedHypoElastic.initializeYourself()

    def setCharacteristicElementLength(self, double length):
        pass

    def computeYourself(self,
                         double[::1] Ke,
                         double[::1] Pe,
                         const double[::1] QTotal,
                         const double[::1] dQ,
                         const double[::1] time,
                         double dTime):

        cdef double pNewDT
        pNewDT = 1e36

        cdef double[::1] S =np.zeros(6)
        cdef double[36] C
        cdef double[6] dKLocal_dDE = np.zeros(6)
        cdef double[6] dS_dK  = np.zeros(6)
        cdef const double[::1] dStrain = dQ[0:6]
        cdef double dK = dQ[self.nU-1]
        cdef double KOld = QTotal[self.nU-1] - dK

        # dummy value for nonLocalRadius
        cdef double nonLocalRadius = 1.0

        S[:] = self.stressLikeInStateVars[0:6]
        cdef double KLocal = self.stressLikeInStateVars[self.nU-1]

        self.marmotMaterialGradientEnhancedHypoElastic.computeStress(
                &S[0],
                KLocal,
                nonLocalRadius,
                &C[0],
                &dKLocal_dDE[0],
                &dS_dK[0],
                &dStrain[0],
                KOld,
                dK,
                &time[0],
                dTime,
                pNewDT)

        Pe[0:6] = S
        Pe[6] = - KLocal + QTotal[6]

        cdef int i = 0
        for i in range(self.nU):
            self.strainLikeInStateVars[i] += dQ [i]

        self.stressLikeInStateVars[0:6] = S
        self.stressLikeInStateVars[self.nU-1] = KLocal

        Ke_ = np.zeros((7,7))

        Ke_[0:6,0:6] = np.array( C ).reshape( (6,6) )
        Ke_[0:6,6] = np.array( dS_dK ).reshape( 6 )
        Ke_[6,0:6] = np.array( dKLocal_dDE ).reshape( 6 )
        Ke_[6,6] = 1.0

        Ke_ravel = Ke_.ravel()
        for i in range(self.nU*self.nU):
            Ke[i] = Ke_ravel[i]
            self.algorithmicTangentInStateVars[i] = Ke_ravel[i]


        if pNewDT < 1.0:
            raise CutbackRequest("Material requests for a cutback!", pNewDT)

    def getNumberOfRequiredStateVars(self,):

        numberOfRequiredStateVarsMaterial = self.marmotMaterialGradientEnhancedHypoElastic.getNumberOfRequiredStateVars()
        numberOfRequiredStateVarsOverhead = 2 * self.nU + self.nU * self.nU
        return numberOfRequiredStateVarsOverhead + numberOfRequiredStateVarsMaterial

    def assignStateVars(self, stateVars):

        self.stateVars = stateVars
        self.stressLikeInStateVars = self.stateVars[0:self.nU]
        self.strainLikeInStateVars = self.stateVars[self.nU:2*self.nU]
        self.algorithmicTangentInStateVars = self.stateVars[2*self.nU:2*self.nU+self.nU*self.nU ]
        self.stateVarsMaterial = self.stateVars[2*self.nU+self.nU*self.nU:]

        self.marmotMaterialGradientEnhancedHypoElastic.assignStateVars(&self.stateVarsMaterial[0],
                                                       len(self.stateVarsMaterial))

    def getResultArray(self, result, getPersistentView=True):

        if result == "stress":
            return np.array(self.stressLikeInStateVars[0:6], copy= not getPersistentView)

        if result == "strain":
            return np.array(self.strainLikeInStateVars[0:6], copy= not getPersistentView)

        if result == "algorithmicTangent":
            return np.array(self.algorithmicTangentInStateVars, copy= not getPersistentView)

        cdef string result_ =  result.encode('UTF-8')

        cdef StateView res = self.marmotMaterialGradientEnhancedHypoElastic.getStateView(result_)

        cdef double[::1] theView = <double[:res.stateSize]> ( res.stateLocation )

        return np.array(  theView, copy= not getPersistentView)

    def __dealloc__(self):

        del self.marmotMaterialGradientEnhancedHypoElastic
