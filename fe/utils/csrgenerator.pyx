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
        
        self.x  = np.zeros_like ( V , dtype=long )
        
        self.data = self.csrMatrix.data
        
        cdef long idx = 0
        cdef long tmp = 0
        cdef long row, col
        cdef int[::1] colIndices, indices, indptr
        indices = self.csrMatrix.indices
        indptr = self.csrMatrix.indptr
        
        for (row, col) in zip(I, J):
            colIndices = indices[indptr[row] : indptr[row+1]]
            tmp = 0 
            for colIdx in colIndices:
                if colIdx == col:
                    self.x [ idx ] = indptr[row] + tmp 
                    break
                tmp += 1
            idx += 1
            
    def updateCSR(self, double[::1] V):
        cdef long i, x
        self.data[:] = 0.0
        for i, x in enumerate( self.x ):
            self.data [ x ] += V [ i ]
        
        return self.csrMatrix.copy()
