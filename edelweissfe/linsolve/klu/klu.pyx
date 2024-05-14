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
#  Alexander Dummer alexander.dummer@uibk.ac.at
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
This module provides an interface to the SuiteSparse KLU solver.
"""

import numpy as np

cimport numpy as np


cdef extern from "kluInterface.h" nogil:
    int solve(double* LHS, int* innerIndices,
              int* outerIndexPtr, int n, double* rhs, int nRhs)


def kluSolve(A, b):
    """
    Solve the linear system A*x = b using the SuiteSparse KLU solver.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        The system matrix.

    b : numpy.ndarray
        The right hand side.

    Returns
    -------
    numpy.ndarray
        The solution x.
    """

    # prepare system matrix
    A = A.tocsc()  # workaround because A is in CSR format and we need CSC

    cdef int rows = A.shape[0]
    cdef double[::1] data = A.data
    cdef int[::1] indices = A.indices
    cdef int[::1] indptr = A.indptr

    # prepare right hand side
    cdef double[::1, :] x = b.reshape((rows, -1), order="F")
    cdef int nRhs = x.shape[1]

    # solution x is stored in b (overwriting b)
    solve(&data[0], &indices[0], &indptr[0], rows, &x[0, 0], nRhs)

    return np.reshape(x, b.shape)
