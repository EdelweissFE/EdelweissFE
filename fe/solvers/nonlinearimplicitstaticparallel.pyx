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
        
        cdef int i, elSizeKe, elNDofPerEl, elNumber, elIdxInVIJ, threadID      
        cdef int desiredThreads = self.numThreads
        cdef int nElements = len(self.elements.values())
        cdef list elList = list(self.elements.values())
        
        cdef double[:, ::1] pNewDTVector = np.ones( (desiredThreads, 1), order='C' )  * 1e36 # as many pNewDTs as threads
        cdef double[:, ::1] Pe = np.empty((desiredThreads,  self.maximumNDofPerEl), ) # Table of Pes - one per thread. 
        cdef double[:, ::1] UN1_ = np.empty((desiredThreads,  self.maximumNDofPerEl), ) # Table of Pes - one per thread. 
        cdef double[:, ::1] dU_ = np.empty((desiredThreads,  self.maximumNDofPerEl), ) # Table of Pes - one per thread. 
        # its size is determined by the 'largest' element (nDofPerEL)
        # note that no Pe can be created inside prange - thus it has to be prepared here
        
        for i in prange(nElements, 
                        schedule='static', 
                        num_threads=desiredThreads, 
                        nogil=True):
            
            threadID = threadid()
            
            with gil:
                el = elList[i]                                                  # python list access - slow
                elIdxInVIJ = self.elementToIndexInVIJMap[el]                    # python dict access - slow
                elSizeKe = el.sizeKe
                elNDofPerEl = el.nDofPerEl
                Pe[threadID, :] = 0.0                                           # prep with zero .. not really necessary
                   
                # we expect that the element to be called releases the GIL!
                el.computeYourself(
                        V   [elIdxInVIJ : elIdxInVIJ + elSizeKe],               # memory view access - fast
                        Pe  [threadID, 0 : elNDofPerEl],                        # memory view access - fast
                        UN1 [ I[ elIdxInVIJ : elIdxInVIJ + elNDofPerEl ]],      # adv. numpy  access - slow
                        dU  [ I[ elIdxInVIJ : elIdxInVIJ + elNDofPerEl ]],      # adv. numpy  access - slow
                        time, dT, 
                        pNewDTVector[threadID, :])                              # memory view access - fast
                
                if pNewDTVector[threadID, 0] <= 1.0:
                    break 
                
                # global effort vector is assembled directly
                P[ I[elIdxInVIJ : elIdxInVIJ+elNDofPerEl] ] += Pe[threadID, 0 : elNDofPerEl]
        
        pNewDT[0] = np.min(pNewDTVector)      

        return P, V, pNewDT