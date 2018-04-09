#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 10:44:29 2018

@author: matthias
"""

from scipy.sparse import csr_matrix
import numpy as np
cimport numpy as np

cdef class CSRGenerator:
    
    cdef csrMatrix
    cdef double[::1] V
    cdef long[::1] x
    cdef long[::1] I, J
    
    cdef double[::1] data
    cdef int nDof
    
    def __init__(self, double[::1] V, long[::1] I, long[::1] J, long nDof):
        
        self.nDof = nDof
        self.I = I
        self.J = J
        self.csrMatrix = csr_matrix( (V, (I,J) ), (nDof, nDof) )        
        self.x  = np.zeros_like ( V , dtype=int )

        cdef int[::1] indices, indptr
        
        self.data = self.csrMatrix.data
        indices = self.csrMatrix.indices
        indptr = self.csrMatrix.indptr
        
        cdef int cooPairIdx, c
        cdef long row, col
        for cooPairIdx in range(len(I)):
            row = I [ cooPairIdx ]
            col = J [ cooPairIdx ]
            for c in range( indptr[row], indptr[row+1] ):
                colIdx = indices[ c ]
                if colIdx == col:
                    self.x [ cooPairIdx ] = c 
                    break
            
    def updateCSR(self, double[::1] V):
        cdef int cooPairIdx
        self.data[:] = 0.0
        for cooPairIdx in range( len (self.x) ):
            self.data [ self.x[cooPairIdx] ] += V [ cooPairIdx ]
        
        return self.csrMatrix.copy()
