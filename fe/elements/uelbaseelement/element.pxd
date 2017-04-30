#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

from fe.materials.umatlibrary cimport pUmatType
from fe.config.ueltypedefs cimport pSimpleUelWithUmatType

cdef extern from "UelInterfaceElement.h":
    cdef cppclass UelInterfaceElement:
        UelInterfaceElement(int elNumber, 
                                        const double* coordinates, 
                                        double *stateVars,
                                        const int nStateVars, 
                                        const double* properties,
                                        int nProperties,
                                        const int* intProperties,
                                        int nIntProperties,
                                        pUmatType umat,
                                        int nStateVarsUmat,
                                        pSimpleUelWithUmatType uel) nogil

        void computeYourself(double* Pe, double* Ke,  const double* UNew, const double* dU,  const double time[], double dTime, double &pNewDT )
        void acceptLastState()

cdef class BaseElement:
    
    cdef UelInterfaceElement* cppBackendElement
    cdef public nodes, 
    cdef public int elNumber
    
    cdef double[::1] uelProperties, stateVars, stateVarsTemp, nodeCoordinates
    cdef pUmatType umat
    cdef pSimpleUelWithUmatType uel
    cdef int nStateVars, nStateVarsUmat
    cdef int[::1] intProperties
    cdef int numGaussPts, uelID
