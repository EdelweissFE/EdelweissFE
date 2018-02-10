#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 20:38:16 2018

@author: matthias
"""

import numpy as np
cimport numpy as np

cdef extern from 'pardisoInterface.h' namespace 'PardisoInterface' nogil:
    int solve(const double* values, const int* innerIndices, const int* outerIndexPtr, int nnz, int cols, int rows,  const double* rhs, double* out)

def pardisoSolve (A, double[::1] b ):
    cdef int rows, cols
    cdef int nnz
    cdef double[::1] data
    cdef int[::1] indices, indptr
    
    rows, cols = A.shape
    nnz = A.nnz
    data = A.data
    indices = A.indices
    indptr = A.indptr
    
    cdef double[::1] x = np.zeros( (rows) )
    
    solve(&data[0], &indices[0], &indptr[0], nnz, cols, rows, &b[0], &x[0])
    
    return np.asarray(x)
