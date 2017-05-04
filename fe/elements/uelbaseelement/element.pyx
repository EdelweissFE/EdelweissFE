#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

import numpy as np
from fe.materials.umatlibrary cimport pUmatType, getUmat
from fe.config.ueltypedefs cimport pSimpleUelWithUmatType
cimport numpy as np
from libcpp.string cimport string
cimport fe.elements.uelbaseelement.element 

cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return False

cdef extern from "userLibrary.h" namespace "userLibrary" nogil:
    BftBaseUel* createBftBaseUel(int id, const double* elementCoordinates,
                                      double* stateVars,
                                              int nStateVars,
                                            const double* properties,
                                            int nProperties,
                                            int noEl,
                                            const pUmatType umat,
                                            int nStateVarsUmat)

cdef class BaseElement:
    
    def __init__(self, nodes, elNumber, int nGaussPt, int uelID):
        self.nodes = nodes
        self.nodeCoordinates = np.concatenate([ node.coordinates for node in nodes])
        self.elNumber = elNumber
        self.numGaussPts = nGaussPt
        self.uelID = uelID
        
    def setProperties(self, uelProperties, umatName, nStateVarsUmat):
        self.uelProperties = uelProperties
        self.nStateVarsUmat = nStateVarsUmat
        self.nStateVars = self.numGaussPts * (nStateVarsUmat + 12)
        self.stateVars = np.zeros(self.nStateVars)
        self.stateVarsTemp = np.zeros(self.nStateVars)
        self.umat = getUmat(umatName.lower())
        self.intProperties = np.empty(0, dtype=np.intc)
        
        if self.bftBaseUel != NULL:
            del self.bftBaseUel
        
        self.bftBaseUel = createBftBaseUel(self.uelID, 
                                                &self.nodeCoordinates[0], 
                                                &self.stateVarsTemp[0], 
                                                self.nStateVars,
                                                &self.uelProperties[0], 
                                                self.uelProperties.shape[0],
                                                self.elNumber,
                                                self.umat, 
                                                self.nStateVarsUmat,)
        if self.bftBaseUel == NULL:
            print("Element not found: {:}".format(self.uelID))

    def initializeStateVarsTemp(self, ):
        self.stateVarsTemp[:] = self.stateVars

    def computeYourself(self, 
                         double[::1] Ke, 
                         double[::1] Pe, 
                         const double[::1] U, 
                         const double[::1] dU, 
                         const double[::1] time, 
                         double dTime, 
                         double[::1] pNewdT):
        self.initializeStateVarsTemp()
        self.bftBaseUel.computeYourself(&U[0], &dU[0],
                                            &Pe[0], &Ke[0],
                                             &time[0],
                                            dTime,  pNewdT[0])
    def acceptLastState(self,):
        self.stateVars[:] = self.stateVarsTemp
        
    def resetToLastValidState(self,):
        pass
    
    resultIndices = {'stress': lambda nStateVarsUmat,kw : slice( kw['gaussPt']*nStateVarsUmat,
                                                                 kw['gaussPt']*nStateVarsUmat + 6) ,
                    'strain': lambda nStateVarsUmat,kw, : slice( kw['gaussPt']*nStateVarsUmat + 6,
                                                                 kw['gaussPt']*nStateVarsUmat + 12),
                    'sdv':    lambda nStateVarsUmat,kw, : slice( kw['idxStart'] + (kw['gaussPt']-1)* nStateVarsUmat, 
                                                                 kw['idxStop' ] + (kw['gaussPt']-1)* nStateVarsUmat)}
    def getPermanentResultPtr(self, **kw):    
        sVars = np.asarray(self.stateVars)
        idxSlice = self.resultIndices[ kw['result'] ](self.nStateVarsUmat, kw)     
        return sVars[idxSlice]
    
    def __dealloc__(self):
        del self.bftBaseUel
