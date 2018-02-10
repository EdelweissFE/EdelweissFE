#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 20:38:16 2018

@author: matthias
"""

import numpy as np
cimport numpy as np

cdef extern from 'pardisoInterface.h' namespace 'PardisoInterface' nogil:
    int solve(const double* values, const int* innerIndices, const int* outerIndexPtr, int nnz, int cols, int rows,  const double* rhs, int nRhs, double* out)

def pardisoSolve (A, b ):
    #A
    cdef int rows, cols 
    rows, cols =            A.shape
    cdef int nnz =          A.nnz
    cdef double[::1] data = A.data
    cdef int[::1] indices = A.indices
    cdef int[::1] indptr =  A.indptr
    
    #b
    cdef double[::1, :] b_ = b.reshape( (rows, -1), order='F') # ensure that we are working on matrices rather than vectors
    cdef int nRhs =         b_.shape[1]
    
    #x
    cdef double[::1, :] x = np.zeros_like(b_, order='F')
    
    solve(&data[0], &indices[0], &indptr[0], nnz, cols, rows, &b_[0,0], nRhs, &x[0,0])
    
    return np.reshape( x, b.shape ) # pass back matrix OR vector, dependent on input
