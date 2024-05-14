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
This module provides an interface to the commercial Panua PARDISO solver.
It is only available, if the PARDISO solver is installed on the system.
You can get the binaries from https://panua.ch/. A license is required to use the solver.
"""

import os

import numpy as np

cimport numpy as np

from edelweissfe.linsolve.pardiso.pardiso import setParameter


cdef extern from "pardiso.h" nogil:
    void pardiso(void*, int*, int *, int*, int *, int *,
                 void*, int*, int *, int*, int *, int *,
                 int*, void*, void*, int*, double*)

    void pardisoinit(void*, int*, int*, int*, double*, int*)


def panuaPardisoSolve(A, b):
    """
    Solve the linear system Ax = b using the Panua PARDISO solver.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        The matrix A of the linear system.
    b : numpy.ndarray
        The right-hand side of the linear system.

    Returns
    -------
    numpy.ndarray
        The solution x of the linear system.
    """

    # prepare matrix
    cdef int rows = A.shape[0]
    cdef double[::1] data = A.data
    cdef int[::1] indices = A.indices + 1  # pardiso uses fortran 1-based indexing
    cdef int[::1] indptr = A.indptr + 1  # pardiso uses fortran 1-based indexing

    # prepare rhs
    cdef double[::1, :] b_ = b.reshape((rows, -1), order="F")
    cdef int nRhs = b_.shape[1]

    # prepare solution vector
    cdef double[::1, :] x = np.zeros_like(b_, order="F")

    # initialize solver
    cdef long int *pt[64]  # internal solver memory pointer
    cdef int[64] iparm  # parameters for pardiso
    cdef double[64] dparm  # parameters for pardiso
    cdef int mtype = 11  # real and unsymmetric matrix
    cdef int solver = 0  # use sparse direct solver
    cdef int error = 0  # initialize error flag

    # check license and initialize the solver
    os.environ["PARDISOLICMESSAGE"] = "1"  # don't print license message
    pardisoinit(pt, &mtype, &solver, iparm, dparm, &error)

    # set parameters
    cdef int maxfct = 5  # maximum number of numerical factorizations
    cdef int mnum = 1  # which factorization to use
    cdef int msglvl = 0  # print statistical information
    cdef int idum = 0  # integer dummy variable
    cdef double ddum = 0  # double dummy variable

    # set custom parameters for pardiso
    iparm = setParameter(iparm, 0, 0)  # use default values

    omp_num_threads = os.environ.get("OMP_NUM_THREADS")
    threads = int(omp_num_threads) if omp_num_threads is not None else -1

    iparm = setParameter(iparm, 2, threads)  # number of threads

    cdef int phase
    # reordering and symbolic factorization
    phase = 11
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &rows, &data[0], &indptr[0], &indices[0], &idum, &nRhs,
            iparm, &msglvl, &ddum, &ddum, &error, dparm)

    # numerical factorization
    phase = 22
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &rows, &data[0], &indptr[0], &indices[0], &idum, &nRhs,
            iparm, &msglvl, &ddum, &ddum, &error, dparm)

    # back substitution and iterative refinement
    phase = 33
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &rows, &data[0], &indptr[0], &indices[0], &idum, &nRhs,
            iparm, &msglvl, &b_[0, 0], &x[0, 0], &error, dparm)

    # free the memory
    phase = -1
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &rows, &data[0], &indptr[0], &indices[0], &idum, &nRhs,
            iparm, &msglvl, &b_[0, 0], &x[0, 0], &error, dparm)

    return np.reshape(x, b.shape)
