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
This module provides an interface to the MUMPS solver.
"""
import numpy as np
from mumps import DMumpsContext


def mumpsSolve(A, b):
    """
    Solve the linear system Ax = b using the MUMPS solver.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        The matrix A.
    b : numpy.ndarray
        The right-hand side vector b.

    Returns
    -------
    numpy.ndarray
        The solution vector x.
    """

    ctx = DMumpsContext()
    ctx.set_silent()

    # analysis step
    if ctx.myid == 0:
        ctx.set_centralized_sparse(A)  # 1-based indexing
    ctx.run(job=1)

    # factorization step
    if ctx.myid == 0:
        ctx.set_centralized_assembled_values(A.data)
    ctx.run(job=2)

    x = np.reshape(b.copy(), (A.shape[0], -1), order="F")
    # solve for multiple right-hand sides
    for i in range(x.shape[1]):
        rhs = x[:, i].copy()
        if ctx.myid == 0:
            ctx.set_rhs(rhs)
        ctx.run(job=3)
        x[:, i] = rhs

    ctx.destroy()  # Cleanup

    return np.reshape(x, b.shape)
