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
#  Daniel Reitmair daniel.reitmair@uibk.ac.at
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

import numpy as np
import pyamg as pa
import scipy.sparse as spa

"""
This module provides an interface to the gmres solver.
"""


def convertListsToTuples(obj):
    if isinstance(obj, dict):
        return {k: convertListsToTuples(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return tuple(convertListsToTuples(x) for x in obj)
    else:
        return obj


class Gmres:
    def __init__(self, opts: dict):
        """Initialize gmres with options.

        Parameters
        ----------
        opts
            File containing the dictionary of options.
        """

        if not opts:
            opts = {}
        opts = convertListsToTuples(opts)
        self.precondOpts = (
            opts["precondopts"]
            if "precondopts" in opts
            else {
                "presmoother": ("block_gauss_seidel", {"iterations": 5}),
                "postsmoother": ("block_gauss_seidel", {"iterations": 5}),
            }
        )
        self.linSolveOpts = opts["linsolveopts"] if "linsolveopts" in opts else {"maxiter": 1, "restart": 500}

    def gmresSolve(self, A: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        Solve the linear system Ax = b using gmres.

        Parameters
        ----------
        A
            Matrix describing the equation system.
        b
            The right hand side vector b.

        Returns
        -------
        x: np.ndarray
            The solution vector x.
        """

        ml = pa.smoothed_aggregation_solver(A, **self.precondOpts)  # get preconditioner
        x, exitCode = spa.linalg.gmres(A, b, atol=1e-9, M=ml.aspreconditioner(), **self.linSolveOpts)
        return np.reshape(x, b.shape)
