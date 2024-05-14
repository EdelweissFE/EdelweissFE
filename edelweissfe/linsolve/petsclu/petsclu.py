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
This module provides an interface to the PetSc LU solver provided by petsc4py.
"""
import numpy as np
from petsc4py import PETSc


def petscluSolve(A, b):
    """
    Solve the linear system Ax = b using the PETSc LU solver.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        The matrix A in the linear system Ax = b.
    b : numpy.ndarray
        The right-hand side vector b in the linear system Ax = b.

    Returns
    -------
    numpy.ndarray
        The solution vector x.
    """

    # Create a PETSc matrix from the SciPy matrix
    # think about not constructing the csr-matrix first and convert it back to aij
    matrix = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))

    ksp = PETSc.KSP().create()
    ksp.setOperators(matrix)
    ksp.setType("preonly")
    ksp.setConvergenceHistory()
    ksp.getPC().setType("lu")

    # prepare rhs
    rhs = b.reshape(A.shape[0], -1)
    nRhs = rhs.shape[1]

    # prepare vectors
    x = matrix.createVecLeft()
    currrhs = matrix.createVecRight()

    # prepare solution array for all rhs
    solution = np.zeros_like(rhs)

    # loop over number of rhs
    for i in range(nRhs):
        # set current rhs
        currrhs.setArray(rhs[:, i])
        # solve for cuurent rhs
        ksp.solve(currrhs, x)
        # store solution
        solution[:, i] = x.getArray()

    ksp.destroy()

    return np.reshape(solution, b.shape)
