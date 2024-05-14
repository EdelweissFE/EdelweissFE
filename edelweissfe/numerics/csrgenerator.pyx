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
# Created on Mon Apr  9 10:44:29 2018

# @author: matthias

import numpy as np
from scipy.sparse import csr_matrix

from edelweissfe.numerics.dofmanager import VIJSystemMatrix

cimport numpy as np


cdef class CSRGenerator:

    cdef csrMatrix
    cdef int[:] x
    cdef double[::1] data
    cdef int nDof, nCooPairs

    def __init__(self, systemMatrix: VIJSystemMatrix):
        """This cdef class generates Compressed Sparse Row Matrices from the COO format,
        and offers the possibility to update the matrix without reanalyzing the
        pattern (in contrast to SciPy).
        Very fast and convenient!

        Parameters
        ----------
        systemMatrix
            The system matrix in COO format, for the first initialization V can be a dummy (empty) vector.
        """

        cdef:
            long[::1] I = systemMatrix.I
            long[::1] J = systemMatrix.J
            long nDof = systemMatrix.nDof

        self.nDof = nDof
        self.x  =  np.zeros_like ( I , dtype=np.intc )

        # abuse Scipy to create the object
        self.csrMatrix = csr_matrix( (np.zeros_like(I, dtype=np.double), (I,J) ), (nDof, nDof) )

        if not self.csrMatrix.has_canonical_format:
            self.csrMatrix.sum_duplicates()

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
            delta = max ( (indptr[row+1] - indptr[row]) >> 1 , 1 )
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

    def updateCSR(self, double[::1] V) -> csr_matrix:
        """Get updated copies of the CSR matrix.

        Returns
        -------
        csr_matrix
            The system matrix in CSR format.
        """

        cdef int cooPairIdx
        self.data[:] = 0.0

        for cooPairIdx in range( self.nCooPairs ):
            self.data [ self.x[cooPairIdx] ] += V [ cooPairIdx ]

        return self.csrMatrix.copy()
