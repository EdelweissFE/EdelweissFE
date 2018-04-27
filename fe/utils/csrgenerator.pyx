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
    cdef int[:] x
    cdef double[::1] data
    cdef int nDof, nCooPairs
    
    def __init__(self, long[::1] I, long[::1] J, long nDof):
        """ Initialize the pattern, V can be a dummy (empty) vector """
        self.nDof = nDof
        self.x  =  np.zeros_like ( I , dtype=np.intc ) 
        
        # abuse Scipy to create the object
        self.csrMatrix = csr_matrix( (np.zeros_like(I, dtype=np.double), (I,J) ), (nDof, nDof) )   

        cdef int[:] indices, indptr
        
        self.data = self.csrMatrix.data
        indices = self.csrMatrix.indices
        indptr = self.csrMatrix.indptr
        
        cdef int cooPairIdx, c, delta
        cdef long row, col
        self.nCooPairs = len(I)
        for cooPairIdx in range(self.nCooPairs):
            row = I [ cooPairIdx ]
            col = J [ cooPairIdx ]
            
            #binary search algorithm (can be improved based on lookup tables, but seems to be sufficient for the moment)
            delta = (indptr[row+1] - indptr[row]) / 2 
            c = indptr[row] + delta
            while True:
                if indices[c] > col:
                    delta = max ( delta >> 1 , 1 )
                    c -= delta
                elif indices[c] < col:
                    delta = max ( delta >> 1 , 1 )
                    c += delta
                else:
                    self.x [ cooPairIdx ] = c 
                    break
            
    def updateCSR(self, double[::1] V):
        """ get updated copies of the CSR matrix """
        cdef int cooPairIdx
        self.data[:] = 0.0
        
        for cooPairIdx in range( self.nCooPairs ):
            self.data [ self.x[cooPairIdx] ] += V [ cooPairIdx ]
            
        return self.csrMatrix.copy()