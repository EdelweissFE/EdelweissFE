#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

from fe.materials.umatlibrary cimport pUmatType
from fe.config.ueltypedefs cimport pSimpleUelWithUmatType

cdef extern from "NISTParallelizableBackendElement.h":
    cdef cppclass NISTParallelizableBackendElement:
        NISTParallelizableBackendElement(int elNumber, 
                                        const double* coordinates, 
                                        double *stateVars,
                                        const int nStateVars, 
#                                        double *stateVarsTemp,
                                        const double* properties,
                                        int nProperties,
                                        const int* intProperties,
                                        int nIntProperties,
                                        pUmatType umat,
                                        int nStateVarsUmat,
                                        pSimpleUelWithUmatType uel) nogil

        void computeYourself(double* Pe, double* Ke,  const double* UNew, const double* dU,  const double time[], double dTime, double &pNewDT )
        void acceptLastState()
        void testIt(int threadId)

cdef class BackendedElement:
    cdef NISTParallelizableBackendElement* backendElement
    cdef NISTParallelizableBackendElement* getBackendElement(self,)
#        return self.