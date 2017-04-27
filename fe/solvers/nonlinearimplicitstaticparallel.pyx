#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 20:37:35 2017

@author: matthias
"""
from fe.solvers.nonlinearimplicitstatic import NIST

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from fe.utils.incrementgenerator import IncrementGenerator
from fe.config.phenomena import flowCorrectionTolerance, effortResidualTolerance
from cython.parallel cimport parallel, threadid, prange
from fe.elements.AbstractBaseElements.CPPBackendedElement cimport BackendedElement
from libc.stdlib cimport malloc, free
from libcpp.string cimport string

cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return False


cdef extern from "NISTParallelizableBackendElement.h":
    cdef cppclass NISTParallelizableBackendElement nogil:
        void computeYourself(double* Pe, double* Ke,  const double* UNew, const double* dU,  const double time[], double dTime, double &pNewDT )

class NISTParallel(NIST):
    """ This is the Nonlinear Implicit STatic -- solver ** Parallel version**.
    Designed to interface with Abaqus UELs
    Public methods are: __init__(), initializeUP() and solveStep(...).
    OutputManagers are updated at the end of each increment. """
    
    identification = "NISTPSolver"
    
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
        
        cdef int elNDofPerEl, elNumber, elIdxInVIJ, threadID      
        cdef int desiredThreads = self.numThreads
        cdef int nElements = len(self.elements.values())
        cdef list elList = list(self.elements.values())
        
        cdef double[:, ::1] pNewDTVector = np.ones( (desiredThreads, 1), order='C' )  * 1e36 # as many pNewDTs as threads
        
        cdef double[::1] UN1_ = UN1
        cdef double[::1] dU_ = dU 
        
        cdef double[:, ::1] UN1e = np.empty((desiredThreads,  self.maximumNDofPerEl), )
        cdef double[:, ::1] dUe = np.empty((desiredThreads,  self.maximumNDofPerEl), )
        cdef double[::1] Pe = np.zeros_like(V) # Table of Pes - one per thread. 
        # its size is determined by the 'largest' element (nDofPerEL)
        # note that no Pe can be created inside prange - thus it has to be prepared here
        
        cdef BackendedElement backendBasedCythonElement
        cdef NISTParallelizableBackendElement** cppBackendElements = <NISTParallelizableBackendElement**> malloc ( nElements * sizeof(NISTParallelizableBackendElement*) )
        cdef int* elIndicesInVIJ = <int*> malloc(nElements * sizeof (int*) )
        cdef int* elNDofs =        <int*> malloc(nElements * sizeof (int*) )
        
        cdef int i, j
        for j in range(nElements):
            backendBasedCythonElement= elList[j]
            cppBackendElements[j] = backendBasedCythonElement.backendElement
            elIndicesInVIJ[j] = self.elementToIndexInVIJMap[backendBasedCythonElement] 
            elNDofs[j] = backendBasedCythonElement.nDofPerEl 
        
        cdef int currentIdxInVIJ
        for i in prange(nElements, 
                    schedule='dynamic', 
                    num_threads=desiredThreads, 
                    nogil=True):
        
            threadID = threadid()
            elIdxInVIJ = elIndicesInVIJ[i]          
            elNDofPerEl = elNDofs[i]# python dict access - slow
            
            for j in range (elNDofPerEl):
                currentIdxInVIJ = I [ elIndicesInVIJ[i]    +  j ]
                UN1e[threadID, j] = UN1_[ currentIdxInVIJ ]
                dUe[threadID, j] = dU_[ currentIdxInVIJ ]
            
            (<NISTParallelizableBackendElement*>  cppBackendElements[i] )[0].computeYourself( &Pe[elIdxInVIJ],
                                    &V[elIdxInVIJ],
                                    &UN1e[threadID, 0],
                                    &dUe[threadID, 0],
                                    &time[0],
                                    dT,
                                    pNewDTVector[threadID, 0])
            
            if pNewDTVector[threadID, 0] <= 1.0:
                break

        pNewDT[0] = np.min(pNewDTVector) 
        
        if pNewDT[0] >= 1.0:
            for i in range(nElements):
                elIdxInVIJ = elIndicesInVIJ[i]
                elNDofPerEl = elNDofs[i]
                for j in range (elNDofPerEl):
                    P[ I[elIdxInVIJ + j ] ] += Pe[ elIdxInVIJ + j ]
                
        free( elIndicesInVIJ )
        free( elNDofs )

        return P, V, pNewDT