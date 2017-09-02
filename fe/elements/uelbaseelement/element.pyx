#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

import numpy as np
#from fe.materials.umatlibrary cimport pUmatType, getUmat
#from fe.config.ueltypedefs cimport pSimpleUelWithUmatType
cimport numpy as np
from libcpp.string cimport string
cimport fe.elements.uelbaseelement.element 

from libc.stdlib cimport malloc, free
cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return False

cdef extern from "userLibrary.h" namespace "userLibrary" nogil:
    BftUel* UelFactory(int id, 
                       const double* elementCoordinates,
                       double* stateVars,
                       int nStateVars,
                       const double* propertiesElement,
                       int nPropertiesElement,
                       int noEl,
                       const string& materialName,
                       int nStateVarsMaterial,
                       const double* propertiesUmat,
                       int nPropertiesUmat)
    
mapLoadTypes={
        'pressure' : DistributedLoadTypes.Pressure
     }

mapStateTypes={
        'geostatic stress' : StateTypes.GeostaticStress
     }

cdef class BaseElement:
    
    def __init__(self, nodes, elNumber, int nGaussPt, int nStateVarsGaussPt, int uelID):
        self.nodes = nodes
        self.nodeCoordinates = np.concatenate([ node.coordinates for node in nodes])
        self.elNumber = elNumber
        self.numGaussPts = nGaussPt
        self.nStateVarsGaussPt = nStateVarsGaussPt
        self.uelID = uelID
        
    def setProperties(self, elementProperties, materialName, nStateVarsMaterial, materialProperties):
        self.elementProperties =    elementProperties
        self.nStateVarsMaterial =   nStateVarsMaterial
        self.materialProperties =   materialProperties
        self.nStateVars =           self.numGaussPts * (nStateVarsMaterial + self.nStateVarsGaussPt)
        self.stateVars =            np.zeros(self.nStateVars)
        self.stateVarsTemp =        np.zeros(self.nStateVars)
        self.materialName =         materialName.upper().encode('UTF-8')
        
        # if we store already an element, we delete it
        if self.bftUel != NULL:
            del self.bftUel
        
        self.bftUel = UelFactory(self.uelID, 
                                                &self.nodeCoordinates[0], 
                                                &self.stateVarsTemp[0], 
                                                self.nStateVars,
                                                &self.elementProperties[0], 
                                                self.elementProperties.shape[0],
                                                self.elNumber,
                                                self.materialName, 
                                                self.nStateVarsMaterial,
                                                &self.materialProperties[0],
                                                self.materialProperties.shape[0])
        if self.bftUel == NULL:
            raise Exception("Element not found: {:}".format(self.uelID))

    def initializeStateVarsTemp(self, ):
        self.stateVarsTemp[:] = self.stateVars
        
    def setInitialCondition(self, 
                            stateType,
                            const double[::1] values):
        
        self.initializeStateVarsTemp()
        self.bftUel.setInitialConditions(mapStateTypes[stateType], &values[0])
        self.acceptLastState()

    def computeYourself(self, 
                         double[::1] Ke, 
                         double[::1] Pe, 
                         const double[::1] U, 
                         const double[::1] dU, 
                         const double[::1] time, 
                         double dTime, 
                         double[::1] pNewdT):
        
        self.initializeStateVarsTemp()
        self.bftUel.computeYourself(&U[0], &dU[0],
                                            &Pe[0], 
                                            &Ke[0],
                                            &time[0],
                                            dTime,  
                                            pNewdT[0])
        
    def computeDistributedLoad(self,
                               str loadType,
                               double[::1] P,
                               int faceID,
                               const double[::1] load,
                               const double[::1] time,
                               double dTime):
        
        self.bftUel.computeDistributedLoad(mapLoadTypes[loadType],
                                    &P[0], 
                                    faceID,
                                    &load[0],
                                    &time[0],
                                    dTime)
    def acceptLastState(self,):
        self.stateVars[:] = self.stateVarsTemp
        
    def resetToLastValidState(self,):
        pass
    
    def getPermanentResultPtr(self, **kw):    
        cdef int gPt =          kw['gaussPt']
        cdef string result =    kw['result'].encode('UTF-8')
        cdef int resultLength = 0
        cdef double* ptr =      self.bftUel.getPermanentResultPointer(result, gPt, resultLength)
        
        if resultLength <= 0:
            raise Exception("Invalid result '{:}' requested for element {:}".format(kw['result'], self.elNumber))
        
        resultArray = np.asarray( <double [:resultLength]> ptr )
        return resultArray[ kw['index'] ] if 'index' in kw else resultArray
    
    def __dealloc__(self):
        del self.bftUel
