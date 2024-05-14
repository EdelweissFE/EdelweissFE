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
Created on Sun Jan  8 20:37:35 2017

@author: matthias

Fri Oct 6 2018:

    This solver is now deprecated since it is replaced by its sucessor mk2.
    This version (mk1) directly accesses the underlying MarmotElement, and thus
    it is able to completely release the GIL throughout the complete prange loop.
    Naturally, it is compatible only with cython elements based on MarmotElements.
    However, it seems that there is no measurable performance advantage over mk2, which
    is not dependent on a  underlying MarmotElement (and is thus also compatible with python or
    cython elements).
    This solver remains for teaching purposes to demonstrate how to interact with cpp objects.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from edelweissfe.solvers.nonlinearimplicitstaticparallelmk2 import NISTParallel
from edelweissfe.utils.exceptions import CutbackRequest

from cython.parallel cimport parallel, prange, threadid
from libc.stdlib cimport free, malloc
from libcpp.string cimport string

from edelweissfe.elements.marmotelement.element cimport (
    MarmotElement,
    MarmotElementWrapper,
)

import os
from multiprocessing import cpu_count
from time import time as getCurrentTime


class NISTParallelForMarmotElements(NISTParallel):

    identification = "NISTPSolver"

    def computeElements(self, elements, Un1, dU, P, K, F, timeStep):
        tic = getCurrentTime()
        cdef:
            double[::1] time = np.array([timeStep.stepTime, timeStep.totalTime])
            double dT = timeStep.timeIncrement

        cdef:
            int elNDofPerEl, elNumber, elIdxInVIJ, elIdxInPe, threadID, currentIdxInU
            int desiredThreads = self.numThreads
            int nElements = len(elements.values())
            list elList = list(elements.values())

            long[::1] I             = K.I
            double[::1] K_mView     = K
            double[::1] UN1_mView   = Un1
            double[::1] dU_mView    = dU
            double[::1] P_mView     = P
            double[::1] F_mView     = F

            double[:, ::1] pNewDTVector = np.ones( (desiredThreads, 1), order='C' )  * 1e36 # as many pNewDTs as threads

            # oversized Buffers for parallel computing:
            # tables [nThreads x max(elements.ndof) ] for U & dU (can be overwritten during parallel computing)
            maxNDofOfAnyEl      = self.theDofManager.largestNumberOfElNDof
            double[:, ::1] UN1e = np.empty((desiredThreads, maxNDofOfAnyEl), )
            double[:, ::1] dUe  = np.empty((desiredThreads, maxNDofOfAnyEl), )
            # oversized buffer for Pe ( size = sum(elements.ndof) )
            double[::1] Pe = np.zeros(self.theDofManager.accumulatedElementNDof)


            MarmotElementWrapper backendBasedCythonElement
            # lists (cpp elements + indices and nDofs), which can be accessed parallely
            MarmotElement** cppElements =      <MarmotElement**> malloc ( nElements * sizeof(MarmotElement*) )
            int[::1] elIndicesInVIJ         = np.empty( (nElements,), dtype=np.intc )
            int[::1] elIndexInPe            = np.empty( (nElements,), dtype=np.intc )
            int[::1] elNDofs                = np.empty( (nElements,), dtype=np.intc )

            int i,j=0

        for i in range(nElements):
            # prepare all lists for upcoming parallel element computing
            backendBasedCythonElement   = elList[i]
            backendBasedCythonElement._initializeStateVarsTemp()
            cppElements[i]              = backendBasedCythonElement.marmotElement
            elIndicesInVIJ[i]           = self.theDofManager.idcsOfHigherOrderEntitiesInVIJ[backendBasedCythonElement]
            elNDofs[i]                  = backendBasedCythonElement.nDof
            # each element gets its place in the Pe buffer
            elIndexInPe[i] = j
            j += elNDofs[i]

        try:
            for i in prange(nElements,
                        schedule='dynamic',
                        num_threads=desiredThreads,
                        nogil=True):

                threadID =      threadid()
                elIdxInVIJ =    elIndicesInVIJ[i]
                elIdxInPe =     elIndexInPe[i]
                elNDofPerEl =   elNDofs[i]

                for j in range (elNDofPerEl):
                    # copy global U & dU to buffer memories for element eval.
                    currentIdxInU =     I [ elIndicesInVIJ[i] +  j ]
                    UN1e[threadID, j] = UN1_mView[ currentIdxInU ]
                    dUe[threadID, j] =  dU_mView[ currentIdxInU ]

                (<MarmotElement*>
                     cppElements[i] )[0].computeYourself(&UN1e[threadID, 0],
                                                        &dUe[threadID, 0],
                                                        &Pe[elIdxInPe],
                                                        &K_mView[elIdxInVIJ],
                                                        &time[0],
                                                        dT,
                                                        pNewDTVector[threadID, 0])

                if pNewDTVector[threadID, 0] <= 1.0:
                    break

            minPNewDT = np.min(pNewDTVector)
            if minPNewDT < 1.0:
                raise CutbackRequest("An element requests for a cutback", minPNewDT)

            #successful elements evaluation: condense oversize Pe buffer -> P
            P_mView[:] = 0.0
            F_mView[:] = 0.0
            for i in range(nElements):
                elIdxInVIJ =    elIndicesInVIJ[i]
                elIdxInPe =     elIndexInPe[i]
                elNDofPerEl =   elNDofs[i]
                for j in range (elNDofPerEl):
                    P_mView[ I[elIdxInVIJ + j] ] +=      Pe[ elIdxInPe + j ]
                    F_mView[ I[elIdxInVIJ + j] ] += abs( Pe[ elIdxInPe + j ] )
        finally:
            free( cppElements )

        toc = getCurrentTime()
        self.computationTimes['elements'] += toc - tic
        return P, K, F
