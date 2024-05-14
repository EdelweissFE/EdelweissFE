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
# Created on Mon Sep 24 13:52:01 2018

# @author: matthias
"""
Parallel implementation of the NIST solver.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from edelweissfe.solvers.nonlinearimplicitstatic import NIST
from edelweissfe.utils.exceptions import CutbackRequest

from cython.parallel cimport parallel, prange, threadid
from libc.stdlib cimport free, malloc
from libcpp.string cimport string

import os
from multiprocessing import cpu_count
from time import time as getCurrentTime


class NISTParallel(NIST):
    identification = "NISTPSolver"

    def solveStep(self, step, model, fieldOutputController, outputmanagers):
        #determine number of threads
        self.numThreads = cpu_count()

        if 'OMP_NUM_THREADS' in os.environ:
            self.numThreads = int( os.environ ['OMP_NUM_THREADS'] ) # higher priority than .inp settings
        # else:
        #     if "NISTSolver" in step.actions["options"]:
        #         self.numThreads = int(step.actions["options"][self.identification].get('numThreads', self.numThreads))

        self.journal.message('Using {:} threads'.format(self.numThreads), self.identification)
        return super().solveStep(step, model, fieldOutputController, outputmanagers)

    def applyDirichletK(self, K, dirichlets):
        """Apply the dirichlet bcs on the global stiffnes matrix
        Is called by solveStep() before solving the global sys.
        http://stackoverflux.com/questions/12129948/scipy-sparse-set-row-to-zeros

        Cythonized version for speed!

        Parameters
        ----------
        K
            The system matrix.
        dirichlets
            The list of dirichlet boundary conditions.

        Returns
        -------
        VIJSystemMatrix
            The modified system matrix.
        """

        cdef int  i, j
        cdef int [::1] indices_, indptr_,
        cdef long[::1] dirichletIndices
        cdef double[::1] data_

        indices_ = K.indices
        indptr_ = K.indptr
        data_ = K.data

        tic = getCurrentTime()
        for dirichlet in dirichlets: # for each bc
            dirichletIndices = self.findDirichletIndices(dirichlet)
            for i in dirichletIndices: # for each node dof in the BC
                for j in range ( indptr_[i] , indptr_ [i+1] ): # iterate along row
                    if i == indices_ [j]:
                        data_[ j ] = 1.0 # diagonal entry
                    else:
                        data_[ j ] = 0.0 # off diagonal entry

        K.eliminate_zeros()
        toc = getCurrentTime()
        self.computationTimes['dirichlet K'] += toc - tic

        return K

    def computeElements(self, elements, Un1, dU, P, K, F, timeStep):

        tic = getCurrentTime()

        cdef:
            double[::1] time = np.array([timeStep.stepTime, timeStep.totalTime])
            double dT = timeStep.timeIncrement

        cdef:
            int elNDof, elNumber, elIdxInVIJ, elIdxInPe, threadID, currentIdxInU
            int desiredThreads = self.numThreads
            int nElements = len(elements.values())
            list elList = list(elements.values())

            long[::1] I             = K.I
            double[::1] K_mView     = K
            double[::1] UN1_mView   = Un1
            double[::1] dU_mView    = dU
            double[::1] P_mView     = P
            double[::1] F_mView     = F

            # oversized Buffers for parallel computing:
            # tables [nThreads x max(elements.ndof) ] for U & dU (can be overwritten during parallel computing)
            maxNDofOfAnyEl      = self.theDofManager.largestNumberOfElNDof
            double[:, ::1] UN1e = np.empty((desiredThreads, maxNDofOfAnyEl), )
            double[:, ::1] dUe  = np.empty((desiredThreads, maxNDofOfAnyEl), )
            # oversized buffer for Pe ( size = sum(elements.ndof) )
            double[::1] Pe = np.zeros(self.theDofManager.accumulatedElementNDof)

            # lists (indices and nDofs), which can be accessed parallely
            int[::1] elIndicesInVIJ         = np.empty( (nElements,), dtype=np.intc )
            int[::1] elIndexInPe            = np.empty( (nElements,), dtype=np.intc )
            int[::1] elNDofs                = np.empty( (nElements,), dtype=np.intc )

            int i, j = 0

        #TODO: this could be done once (__init__) and stored permanently in a cdef class
        for i in range(nElements):
            # prepare all lists for upcoming parallel element computing
            el                      = elList[i]
            elIndicesInVIJ[i]       = self.theDofManager.idcsOfHigherOrderEntitiesInVIJ[el]
            elNDofs[i]              = el.nDof
            # each element gets its place in the Pe buffer
            elIndexInPe[i] = j
            j += elNDofs[i]

        for i in prange(nElements,
                    schedule='dynamic',
                    num_threads=desiredThreads,
                    nogil=True):

            threadID    = threadid()
            elIdxInVIJ  = elIndicesInVIJ[i]
            elIdxInPe   = elIndexInPe[i]
            elNDof      = elNDofs[i]

            for j in range (elNDof):
                # copy global U & dU to buffer memories for element eval.
                currentIdxInU =     I [ elIndicesInVIJ[i] +  j ]
                UN1e[threadID, j] = UN1_mView[ currentIdxInU ]
                dUe[threadID, j] =  dU_mView[ currentIdxInU ]

            # for accessing the element in the list, and for passing the parameters
            # we have to enable the gil.
            # This prevents a parallel execution in the meanwhile,
            # so we hope the method computeYourself AGAIN releases the gil INSIDE.
            # Otherwise, a truly parallel execution won't happen at all!
            with gil:
                elList[i].computeYourself(K_mView[elIdxInVIJ : elIdxInVIJ + elNDof ** 2],
                                          Pe[elIdxInPe : elIdxInPe + elNDof],
                                          UN1e[threadID, :],
                                          dUe[threadID, :],
                                          time,
                                          dT)

        #successful elements evaluation: condense oversize Pe buffer -> P
        for i in range(nElements):
            elIdxInVIJ =    elIndicesInVIJ[i]
            elIdxInPe =     elIndexInPe[i]
            elNDof =   elNDofs[i]
            for j in range (elNDof):
                P_mView[ I[elIdxInVIJ + j] ] +=      Pe[ elIdxInPe + j ]
                F_mView[ I[elIdxInVIJ + j] ] += abs( Pe[ elIdxInPe + j ] )

        toc = getCurrentTime()
        self.computationTimes['elements'] += toc - tic

        return P, K, F
