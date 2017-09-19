#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 20:57:09 2017

@author: matthias
"""

    
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc

cdef class ElementResultCollector:
    cdef int nGauss, nEls, nSize 
    cdef public  resultsTable
    cdef double[:, :, ::1] res_
    
    cdef double** resultPointers
    
    def __init__(self, elements, gaussPoints, result):
        
        self.nEls = len(elements)
        self.nGauss = len(gaussPoints)
        
        resultsPointerList = [ [ el.getPermanentResultPtr(result = result, gaussPt = gpt) for gpt in gaussPoints ] for el in elements ]
        
        self.nSize = resultsPointerList[0][0].shape[0]
        
        self.resultPointers = <double**> malloc ( sizeof(double*) * self.nEls * self.nGauss )
        
        cdef double* ptr
        cdef double[::1] res
        for i, el in enumerate(resultsPointerList):
            for j, gPt in enumerate(el):
                res = gPt
                ptr = <double*> &res[0]
                self.resultPointers[ i * self.nGauss + j ] = ptr
        
        self.resultsTable = np.empty([self.nEls, self.nGauss, self.nSize]  )
        self.res_ = self.resultsTable
        
        if self.nGauss == 1:
            self.resultsTable = self.resultsTable.reshape(self.nEls, -1)
    
    def update(self, ):
        cdef int i, j, k
        for i in range(self.nEls):
            for j in range(self.nGauss):
                for k in range(self.nSize):
                    self.res_[i,j,k] = self.resultPointers[ i * self.nGauss + j ][k]
    
    def getCurrentResults(self,):
        self.update()
        return self.resultsTable
