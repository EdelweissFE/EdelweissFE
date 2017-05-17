#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

from fe.materials.umatlibrary cimport pUmatType

cdef extern from "bftUel.h" namespace "BftUel":
    cdef enum StateTypes:
        Sigma11,
        Sigma22,
        Sigma33,
        HydrostaticStress,
        GeostaticStress,
        UmatStateVars
        
    cdef enum DistributedLoadTypes:
        Pressure

cdef extern from "bftUel.h":
    cdef cppclass BftUel nogil:
        void computeYourself( const double* QTotal,
                                            const double* dQ,
                                            double* Pe,
                                            double* Ke,
                                            const double* time,
                                            double dT,
                                            double& pNewdT,)
        
        void setInitialConditions(StateTypes state, 
                                  const double* values)
        
        void computeDistributedLoad(
                                DistributedLoadTypes loadType,
                                double* P, 
                                int faceID, 
                                const double* load,
                                const double* time,
                                double dT)
        
cdef class BaseElement:
    
    cdef BftUel* bftUel
    cdef public nodes, 
    cdef public int elNumber
    
    cdef public double[::1] stateVars
    cdef double[::1] uelProperties, stateVarsTemp, nodeCoordinates, umatProperties
    cdef pUmatType umat
    cdef public nStateVars, nStateVarsUmat
    cdef int[::1] intProperties
    cdef int numGaussPts, uelID
