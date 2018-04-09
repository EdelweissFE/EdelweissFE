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
    
    """ Generates Compressed Sparse Row Matrices from the COO format,
    and offers the possibility to update the matrix without reanalyzing the 
    pattern (in contrast to SciPy)"""
    
    cdef csrMatrix
    cdef long[::1] x
    cdef double[::1] data
    cdef int nDof
    
    def __init__(self, long[::1] I, long[::1] J, long nDof):
        """ Initialize the pattern, V can be a dummy (empty) vector """
        self.nDof = nDof     
        self.x  = np.zeros_like ( I , dtype=int )
        
        # abuse Scipy to create the object
        self.csrMatrix = csr_matrix( (np.zeros_like(I, dtype=np.double), (I,J) ), (nDof, nDof) )   

        cdef int[::1] indices, indptr
        
        self.data = self.csrMatrix.data
        indices = self.csrMatrix.indices
        indptr = self.csrMatrix.indptr
        
        cdef int cooPairIdx, c, delta, nCols, base
        cdef long row, col
        
        for cooPairIdx in range(len(I)):
            row = I [ cooPairIdx ]
            col = J [ cooPairIdx ]
            
            #binary search algorithm (can be improved based on lookup table, but seems to be sufficient for the moment)
            base = indptr[row] 
            nCols = indptr[row+1] - indptr[row]
            delta = nCols / 2 
            c = base + delta
            while True:
                colIdx = indices [ c  ]
                if colIdx > col:
                    delta = max ( delta >> 1 , 1 )
                    c -= delta
                elif colIdx < col:
                    delta = max ( delta >> 1 , 1 )
                    c += delta
                else:
                    self.x [ cooPairIdx ] = c 
                    break
            
#            for c in range( indptr[row], indptr[row+1] ):
#                colIdx = indices[ c ]
#                if colIdx == col:
#                    self.x [ cooPairIdx ] = c 
#                    break
                
        
            
    def updateCSR(self, double[::1] V):
        """ get updated copies of the CSR matrix """
        cdef int cooPairIdx
        self.data[:] = 0.0
        for cooPairIdx in range( len (self.x) ):
            self.data [ self.x[cooPairIdx] ] += V [ cooPairIdx ]
        
        return self.csrMatrix.copy()
