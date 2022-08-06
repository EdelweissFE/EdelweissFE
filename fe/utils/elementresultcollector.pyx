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
Collecting elemental results may be a perfomance critical part.
This module provides a cdef class for the efficient gathering.
"""

    
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

cdef class ElementResultCollector:
    """
    A cython class for collecting element results (by using the permanent results pointer (i.e., a numpy array) 
    in large array of all elements and all gaussPoints.
    A 3D array is assembled if multiple gaussPoints are requested (shape [elements, gaussPoints, resultVector] )
    or a 2D array for one gaussPoint ( shape [elements, resultVector] ).
    
    method getCurrentResults() updates the assembly array and passes it back.
    
    The caller is responsible to make a copy of it, if persistent results are needed!
    """
    
    cdef public  resultsTable
    
    cdef int nGauss, nEls, nSize 
    cdef double[:, :, ::1] res_
    cdef double** resultPointers
    
    def __init__(self, elements, gaussPoints, result):
        
        self.nEls = len(elements)
        self.nGauss = len(gaussPoints)
        # assemble a 2d list of all permanent result arrays (=continously updated np arrays)
        resultsPointerList = [ [ el.getResultArray(result, qp, getPersistentView=True) for qp in gaussPoints ] for el in elements ]
        self.nSize = resultsPointerList[0][0].shape[0]
        
        # allocate an equivalent 2D C-array for the pointers to each elements results
        self.resultPointers = <double**> malloc ( sizeof(double*) * self.nEls * self.nGauss )
        
        cdef double* ptr
        cdef double[::1] res
        # fill the 2D C-array of pointers by accessing the elements' resultArrays memoryviews
        for i, el in enumerate(resultsPointerList):
            for j, gPt in enumerate(el):
                res = gPt
                ptr = <double*> &res[0]
                self.resultPointers[ i * self.nGauss + j ] = ptr
        
        # initialize the large assembly array
        self.resultsTable = np.empty([self.nEls, self.nGauss, self.nSize]  )
        
        # an internal use only memoryview is created for accesing the assembly
        self.res_ = self.resultsTable
        
        # remove the useless axis
        if self.nGauss == 1:
            self.resultsTable = self.resultsTable.reshape(self.nEls, -1)
    
    def update(self, ):
        #performant updating!
        cdef int i, j, k
        for i in range(self.nEls):
            for j in range(self.nGauss):
                for k in range(self.nSize):
                    # most inner loop: could also be handled by copying the complete vector at once,
                    # but this version turned out to be faster!
                    self.res_[i,j,k] = self.resultPointers[ i * self.nGauss + j ][k]
    
    def getCurrentResults(self,):
        self.update()
        return self.resultsTable
    
    def __dealloc__(self):
        free ( self.resultPointers )

