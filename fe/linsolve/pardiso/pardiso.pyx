#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  ---------------------------------------------------------------------
#
#  _____    _      _              _         _____ _____ 
# | ____|__| | ___| |_      _____(_)___ ___|  ___| ____|
# |  _| / _` |/ _ \ \ \ /\ / / _ \ / __/ __| |_  |  _|  
# | |__| (_| |  __/ |\ V  V /  __/ \__ \__ \  _| | |___ 
# |_____\__,_|\___|_| \_/\_/ \___|_|___/___/_|   |_____|
#                                                       
# 
#  Unit of Strength of Materials and Structural Analysis
#  University of Innsbruck,
#  2017 - today
# 
#  Matthias Neuner matthias.neuner@uibk.ac.at
# 
#  This file is part of EdelweissFE.
# 
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
# 
#  The full text of the license can be found in the file LICENSE.md at
#  the top level directory of EdelweissFE.
#  ---------------------------------------------------------------------
"""
Created on Fri Feb  9 20:38:16 2018

@author: matthias
"""

import numpy as np
cimport numpy as np

cdef extern from 'pardisoInterface.h' namespace 'PardisoInterface' nogil:
    int solve(const double* values, const int* innerIndices, const int* outerIndexPtr, int nnz, int cols, int rows,  const double* rhs, int nRhs, double* out)

def pardisoSolve (A, b ):
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
