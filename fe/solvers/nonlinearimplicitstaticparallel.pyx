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
from fe.elements.uelbaseelement.element cimport BaseElement, BftUel
from libc.stdlib cimport malloc, free
from libcpp.string cimport string
from time import time as getCurrentTime

cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#     print(cppString.decode('UTF-8'))
    return False

class NISTParallel(NIST):
    """ This is the Nonlinear Implicit STatic -- solver ** Parallel version**.
    Designed to interface with Abaqus UELs
    Public methods are: __init__(), initializeUP() and solveStep(...).
    OutputManagers are updated at the end of each increment. """
    
    identification = "NISTPSolver"
    
    def __init__(self, jobInfo, modelInfo, journal, outputmanagers=None):
        super().__init__(jobInfo, modelInfo, journal, outputmanagers)
        
        self.maximumNDofPerEl = 0
        self.sizePe = 0
        for el in self.elements.values():
            self.maximumNDofPerEl = el.nDofPerEl if el.nDofPerEl > self.maximumNDofPerEl else self.maximumNDofPerEl
            self.sizePe += el.nDofPerEl
            
    def solveStep(self, step, time, stepActions, stepOptions, U, P):
        """ Public interface to solve for an ABAQUS like step
        returns: boolean Success, U vector, P vector, and the new current total time """
        
        self.numThreads = int(stepOptions['NISTSolver'].get('numThreads', 1))
        return super().solveStep(step, time, stepActions, stepOptions, U, P)
    
    def computeElements(self, U, dU, double[::1] time, double dT,
                        pNewDT, 
                        P, 
                        V, 
                        long[::1] I, 
                        long[::1] J):
        """ Loop over all elements, and evalute them. 
        Note that ABAQUS style is employed: element(Un+1, dUn+1) 
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration"""

        UN1 = dU + U # ABAQUS style!
        
        cdef int elNDofPerEl, elNumber, elIdxInVIJ, elIdxInPe, threadID, currentIdxInU   
        cdef int desiredThreads = self.numThreads
        cdef int nElements = len(self.elements.values())
        cdef list elList = list(self.elements.values())
        
        cdef double[::1] V_mView = V
        cdef double[:, ::1] pNewDTVector = np.ones( (desiredThreads, 1), order='C' )  * 1e36 # as many pNewDTs as threads
        cdef double[::1] UN1_mView = UN1
        cdef double[::1] dU_mView = dU 

        # oversized Buffers for parallel computing:
        # tables [nThreads x max(elements.ndof) ] for U & dU (can be overwritten during parallel computing)
        cdef double[:, ::1] UN1e = np.empty((desiredThreads,  self.maximumNDofPerEl), )
        cdef double[:, ::1] dUe = np.empty((desiredThreads,  self.maximumNDofPerEl), )
        # oversized buffer for Pe ( size = sum(elements.ndof) )
        cdef double[::1] Pe = np.empty(self.sizePe) 
        
        cdef BaseElement backendBasedCythonElement
        # lists (cpp elements + indices and nDofs), which can be accessed parallely
        cdef BftUel** cppElements = <BftUel**> malloc ( nElements * sizeof(BftUel*) )
        cdef int* elIndicesInVIJ = <int*> malloc(nElements * sizeof (int*) )
        cdef int* elIndexInPe =    <int*> malloc(nElements * sizeof (int*) )
        cdef int* elNDofs =        <int*> malloc(nElements * sizeof (int*) )
        
        cdef int i, j=0
        for i in range(nElements):
            # prepare all lists for upcoming parallel element computing
            backendBasedCythonElement=  elList[i]
            backendBasedCythonElement.initializeStateVarsTemp()
            cppElements[i] =     backendBasedCythonElement.bftUel
            elIndicesInVIJ[i] =         self.elementToIndexInVIJMap[backendBasedCythonElement] 
            elNDofs[i] =                backendBasedCythonElement.nDofPerEl 
            # each element gets its place in the Pe buffer
            elIndexInPe[i] = j
            j += elNDofs[i]
        
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
            
            (<BftUel*> 
                 cppElements[i] )[0].computeYourself(&UN1e[threadID, 0],
                                                    &dUe[threadID, 0],
                                                    &Pe[elIdxInPe],
                                                    &V_mView[elIdxInVIJ],
                                                    &time[0],
                                                    dT,
                                                    pNewDTVector[threadID, 0])
            
            if pNewDTVector[threadID, 0] <= 1.0:
                break

        pNewDT[0] = np.min(pNewDTVector) 
        
        cdef double[::1] P_mView = P
        if pNewDT[0] >= 1.0:
            #successful elements evaluation: condense oversize Pe buffer -> P
            P_mView[:] = 0.0
            for i in range(nElements):
                elIdxInVIJ =    elIndicesInVIJ[i]
                elIdxInPe =     elIndexInPe[i]
                elNDofPerEl =   elNDofs[i]
                for j in range (elNDofPerEl): 
                    P_mView[ I[elIdxInVIJ + j] ] += Pe[ elIdxInPe + j ]
                
        free( elIndicesInVIJ )
        free( elNDofs )
        free( elIndexInPe )
        free( cppElements )

        return P, V, pNewDT
