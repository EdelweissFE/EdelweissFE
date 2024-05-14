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
# Created on Fri Feb  9 20:38:16 2018

# @author: matthias

"""
This module provides an interface to the PARDISO solver provided by the Intel Math Kernel Library (MKL).
"""

import os

import numpy as np

cimport numpy as np


cdef extern from "mkl.h" nogil:
    void pardiso(void*, int*,  int*, int*,
                 int*,  int*, void*, int*,
                 int*,  int*,  int*, int*,
                 int*, void*, void*, int*)

    void pardisoinit(void*, int*, int*)


def setParameter(iparm, idx: int, value: int):
    """
    Set a parameter of the PARDISO solver.

    Parameters
    ----------
    iparm : int[:]
        The parameter array.
    idx : int
        The index of the parameter.
    value : int
        The value of the parameter.
    """
    iparm[idx] = value

    return iparm


def pardisoSolve(A, b):
    """
    Solve a linear system of equations using the Intel MKL PARDISO solver.

    Parameters
    ----------
    A : csr_matrix
        The system matrix.
    b : ndarray
        The right-hand side vector.

    Returns
    -------
    ndarray
        The solution vector.
    """

    # prepare system matrix
    cdef int rows = A.shape[0]
    cdef double[::1] data = A.data
    cdef int[::1] indices = A.indices + 1  # pardiso uses fortran 1-based indexing
    cdef int[::1] indptr = A.indptr + 1    # pardiso uses fortran 1-based indexing

    # prepare rhs
    cdef double[::1, :] b_ = b.reshape((rows, -1), order="F")
    cdef int nRhs = b_.shape[1]

    # prepare solution vector
    cdef double[::1, :] x = np.zeros_like(b_, order="F")

    # initialize solver
    cdef long int *pt[64]  # internal solver memory pointer
    cdef int[64] iparm     # parameters for pardiso
    cdef int mtype = 11    # real and unsymmetric matrix
    cdef int error = 0     # initialize error flag

    # initialie pardiso solver
    pardisoinit(pt, &mtype, &iparm[0])

    # set parameters
    cdef int maxfct = 5   # maximum number of numerical factorizations
    cdef int mnum = 1     # which factorization to use
    cdef int msglvl = 0   # print statistical information
    cdef int[::1] perm    # permutation vector
    cdef double ddum = 0  # dummy variable

    # set custom parameters for pardiso
    iparm = setParameter(iparm, 0, 0)  # use default values

    omp_num_threads = os.environ.get("OMP_NUM_THREADS")
    threads = int(omp_num_threads) if omp_num_threads is not None else -1

    # set parameters for solver
    # usage: setParameter(iparm, idx, value)
    iparm = setParameter(iparm, 2, threads)  # set number of threads

    cdef int phase
    # reordering and symbolic factorization
    phase = 11
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &rows, &data[0], &indptr[0], &indices[0], &perm[0], &nRhs,
            iparm, &msglvl, &ddum, &ddum, &error)

    # numerical factorization
    phase = 22
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &rows, &data[0], &indptr[0], &indices[0], &perm[0], &nRhs,
            &iparm[0], &msglvl, &ddum, &ddum, &error)

    # back substitution and iterative refinement
    phase = 33
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &rows, &data[0], &indptr[0], &indices[0], &perm[0], &nRhs,
            &iparm[0], &msglvl, &b_[0, 0], &x[0, 0], &error)

    # free the memory
    phase = -1
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &rows, &data[0], &indptr[0], &indices[0], &perm[0], &nRhs,
            &iparm[0], &msglvl, &b_[0, 0], &x[0, 0], &error)

    return np.reshape(x, b.shape)
