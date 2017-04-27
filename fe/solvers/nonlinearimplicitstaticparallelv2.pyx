#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 20:37:35 2017

@author: matthias
"""
from fe.solvers.nonlinearimplicitstatic import NIST

import numpy as np
#cimport numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from fe.utils.incrementgenerator import IncrementGenerator
from fe.config.phenomena import flowCorrectionTolerance, effortResidualTolerance
from cython.parallel cimport parallel, threadid, prange

from fe.elements.AbstractBaseElements.CPPBackendedElement cimport BackendedElement
from libc.stdlib cimport malloc, free

cdef extern from "NISTParallelizableBackendElement.h":
    cdef cppclass NISTParallelizableBackendElement nogil:
        void computeYourself(double* Pe, double* Ke,  const double* UNew, const double* dU,  const double time[], double dTime, double &pNewDT )
        void testIt(int tId)
        
#cdef class BackendedElement:
#    cdef NISTParallelizableBackendElement* backendElement

class NISTParallel(NIST):
    """ This is the Nonlinear Implicit STatic -- solver ** Parallel version**.
    Designed to interface with Abaqus UELs
    Public methods are: __init__(), initializeUP() and solveStep(...).
    OutputManagers are updated at the end of each increment. """
    
    identification = "NISTPSolver"
    
#    cdef int test
    
    def __init__(self, jobInfo, modelInfo, journal, outputmanagers=None):
        super().__init__(jobInfo, modelInfo, journal, outputmanagers)
        
        self.maximumNDofPerEl = 0
        for el in self.elements.values():
            self.maximumNDofPerEl = el.nDofPerEl if el.nDofPerEl > self.maximumNDofPerEl else self.maximumNDofPerEl
            
        
        
    def solveStep(self, step, time, stepActions, U, P):
        """ Public interface to solve for an ABAQUS like step
        returns: boolean Success, U vector, P vector, and the new current total time """
        
        self.numThreads = int(stepActions['NISTSolverOptions'].get('numThreads', 1))
        return super().solveStep(step, time, stepActions, U, P)
    
    def computeElements(self, U, dU, double[::1] time, double dT,
                        double[::1] pNewDT, 
                        P, 
                        double[::1] V, 
                        long[::1] I, 
                        long[::1] J):
        """ Loop over all elements, and evalute them. 
        Note that ABAQUS style is employed: element(Un+1, dUn+1) 
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration 
        UNOPTIMIZED PROTOTYPE."""
            
        P[:] = 0.0
        UN1 = dU + U # ABAQUS style!
        
        cdef int i, elSizeKe, elNDofPerEl, elNumber, elIdxInVIJ, threadID      
        cdef int desiredThreads = self.numThreads
        cdef int nElements = len(self.elements.values())
        cdef list elList = list(self.elements.values())
        
        cdef double[:, ::1] pNewDTVector = np.ones( (desiredThreads, 1), order='C' )  * 1e36 # as many pNewDTs as threads
#        cdef double[:, ::1] Pe = np.empty((desiredThreads,  self.maximumNDofPerEl), ) # Table of Pes - one per thread. 
        
        cdef double[::1] UN1_ = UN1
        cdef double[::1] dU_ = dU 
        
        cdef double[:, ::1] UN1e_ = np.empty((desiredThreads,  self.maximumNDofPerEl), ) # Table of Pes - one per thread. 
        cdef double[:, ::1] dUe_ = np.empty((desiredThreads,  self.maximumNDofPerEl), ) # Table of Pes - one per thread. 
        cdef double[::1] Pe_ = np.zeros_like(V) # Table of Pes - one per thread. 
        # its size is determined by the 'largest' element (nDofPerEL)
        # note that no Pe can be created inside prange - thus it has to be prepared here
        
        cdef BackendedElement el_
        cdef NISTParallelizableBackendElement** elements_ = <NISTParallelizableBackendElement**> malloc ( nElements * sizeof(NISTParallelizableBackendElement*) )
        cdef int* elIndicesInVIJ_ = <int*> malloc(nElements * sizeof (int*) )
        cdef int* elNDofs_ =        <int*> malloc(nElements * sizeof (int*) )
        cdef int j_
            
        
        for j_ in range(nElements):
            el_= elList[j_]
            elements_[j_] = el_.backendElement
            elIndicesInVIJ_[j_] = self.elementToIndexInVIJMap[el_] 
            elNDofs_[j_] = el_.nDofPerEl 
        
        
        cdef int parI_
        
        cdef int innerParI, currentIdxInVIJ
        
        cdef double pnewdt
        for parI_ in prange(nElements, 
                    schedule='dynamic', 
                    num_threads=desiredThreads, 
                    nogil=True):
        
            threadID = threadid()
#            pnewdt = 1e36
            
#            if pnewdt >= 1.0:
#                break
            
            elIdxInVIJ = elIndicesInVIJ_[parI_]          
            elNDofPerEl = elNDofs_[parI_]# python dict access - slow
            elSizeKe = elNDofPerEl* elNDofPerEl
            
            for innerParI in range (elNDofPerEl):
                currentIdxInVIJ = I [ elIndicesInVIJ_[parI_]    +  innerParI ]
                UN1e_[threadID, innerParI] = UN1_[ currentIdxInVIJ ]
                dUe_[threadID, innerParI] = dU_[ currentIdxInVIJ ]
            
#            
            #void computeYourself(double* Pe, double* Ke,  const double* UNew, const double* dU,  const double time[], double dTime, double &pNewDT )
            
#           cdef NISTParallelizableBackendElement* parEl_ 
#            parEl_ = elements_[parI_] --> cython compiler bug
            (<NISTParallelizableBackendElement*>elements_[parI_])[0].computeYourself( &Pe_[elIdxInVIJ],
                                    &V[elIdxInVIJ],
                                    &UN1e_[threadID, 0],
                                    &dUe_[threadID, 0],
                                    &time[0],
                                    dT,
                                    pNewDTVector[threadID, 0]
                                   )
            
            if pNewDTVector[threadID, 0] <= 1.0:
                break

        pNewDT[0] = np.min(pNewDTVector) 
        
        if pNewDT[0] >= 1.0:
            for el in elList:
                elIdxInVIJ = self.elementToIndexInVIJMap[el]                    # python dict access - slow
                elNDofPerEl = el.nDofPerEl
                P[ I[elIdxInVIJ : elIdxInVIJ+elNDofPerEl] ] += Pe_[ elIdxInVIJ : elIdxInVIJ + elNDofPerEl ]

        free( elements_ )
        free( elIndicesInVIJ_ )
        free( elNDofs_ )
             

        return P, V, pNewDT