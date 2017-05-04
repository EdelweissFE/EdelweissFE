#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

from fe.materials.umatlibrary cimport pUmatType

cdef extern from "userLibrary.h":
        cdef cppclass BftBaseUel nogil:
            void computeYourself( const double* QTotal,
                                                const double* dQ,
                                                double* Pe,
                                                double* Ke,
                                                const double* time,
                                                double dT,
                                                double& pNewdT,)
            
cdef class BaseElement:
    
    cdef BftBaseUel* bftBaseUel
    cdef public nodes, 
    cdef public int elNumber
    
    cdef double[::1] uelProperties, stateVars, stateVarsTemp, nodeCoordinates
    cdef pUmatType umat
    cdef int nStateVars, nStateVarsUmat
    cdef int[::1] intProperties
    cdef int numGaussPts, uelID
