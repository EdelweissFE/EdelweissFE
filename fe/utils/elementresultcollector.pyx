#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 20:57:09 2017

@author: matthias
"""

    
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

cdef class ElementResultCollector:
    """
    A cython class for collecting element results (by using the permanent results pointer (i.e., a numpy array) 
    into large array of all elements and all gaussPoints.
    A 3D array is assembled if multiple gaussPoints are requested (shape [elements, gaussPoints, resultVector] )
    or a 2D array for one gaussPoint ( shape [elements, resultVector] ).
    
    Method getCurrentResults updates the assembly array and passes it back.
    
    The caller is responsible to make a copy of it, if persistent results are needed!
    """
    
    cdef public  resultsTable
    
    cdef int nGauss, nEls, nSize 
    cdef double[:, :, ::1] res_
    cdef double** resultPointers
    
    def __init__(self, elements, gaussPoints, result):
        
        self.nEls = len(elements)
        self.nGauss = len(gaussPoints)
        # assemble a 2d list of all permanent result pointers 
        resultsPointerList = [ [ el.getPermanentResultPtr(result = result, gaussPt = gpt) for gpt in gaussPoints ] for el in elements ]
        self.nSize = resultsPointerList[0][0].shape[0]
        
        # use a performant 2D C-array for the pointers
        self.resultPointers = <double**> malloc ( sizeof(double*) * self.nEls * self.nGauss )
        
        cdef double* ptr
        cdef double[::1] res
        # fill the 2D C-array of pointers
        for i, el in enumerate(resultsPointerList):
            for j, gPt in enumerate(el):
                res = gPt
                ptr = <double*> &res[0]
                self.resultPointers[ i * self.nGauss + j ] = ptr
        
        # initialize the large assembly array
        self.resultsTable = np.empty([self.nEls, self.nGauss, self.nSize]  )
        
        # an internal use only memoryview is created for accesing the assembly
        self.res_ = self.resultsTable
        
        # remove the useless axis
        if self.nGauss == 1:
            self.resultsTable = self.resultsTable.reshape(self.nEls, -1)
    
    def update(self, ):
        
        #performant updating!
        cdef int i, j, k
        for i in range(self.nEls):
            for j in range(self.nGauss):
                for k in range(self.nSize):
                    # most inner loop: could also be handled by copying the complete vector at once,
                    # but this version turned out to be faster!
                    self.res_[i,j,k] = self.resultPointers[ i * self.nGauss + j ][k]
    
    def getCurrentResults(self,):
        self.update()
        return self.resultsTable
    
    def __dealloc__(self):
        free ( self.resultPointers )

