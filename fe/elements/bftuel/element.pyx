#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

import numpy as np
cimport numpy as np
cimport fe.elements.bftuel.element 
cimport libcpp.cast

from libcpp.memory cimport unique_ptr, allocator, make_unique

from libc.stdlib cimport malloc, free
cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return False
    
mapLoadTypes={
        'pressure' : DistributedLoadTypes.Pressure
     }

mapStateTypes={
        'geostatic stress' : StateTypes.GeostaticStress
     }
    
cdef class BftUelWrapper:
    
    def __cinit__(self, elementType, nodes, elNumber):
        self.nodes = nodes
        if self.nodes:
            self.nodeCoordinates = np.concatenate([ node.coordinates for node in nodes])
            
        self.elNumber = elNumber
        
        self.bftUel = UelFactory(getElementCodeFromName( elementType.upper().encode('utf-8')), 
                        self.elNumber)
        
        self.nNodes                         = self.bftUel.getNNodes()
        
        self.nDofPerEl                      = self.bftUel.getNDofPerElement()
        
        cdef vector[vector[string]] fields  = self.bftUel.getNodeFields()
        self.fields                         = [ [ s.decode('utf-8')  for s in n  ] for n in fields ]
        
        cdef vector[int] permutationPattern = self.bftUel.getDofIndicesPermutationPattern()
        self.dofIndicesPermutation          = np.asarray(permutationPattern)
        
        self.ensightType                    = self.bftUel.getElementShape().decode('utf-8')
        
    def setProperties(self, elementProperties, materialName, materialProperties):
        
        self.elementProperties =    elementProperties
        self.materialProperties =   materialProperties
        
        self.bftUel.assignProperty( 
                ElementProperties(
                        &self.elementProperties[0],
                        self.elementProperties.shape[0] ) )
        
        self.bftUel.assignProperty(
                BftMaterialSection(
                        getMaterialCodeFromName(
                                materialName.upper().encode('UTF-8')), 
                        &self.materialProperties[0],
                        self.materialProperties.shape[0] ) )
        
        self.nStateVars =           self.bftUel.getNumberOfRequiredStateVars()
        
        self.stateVars =            np.zeros(self.nStateVars)
        self.stateVarsTemp =        np.zeros(self.nStateVars)
        
        self.bftUel.assignStateVars(&self.stateVarsTemp[0], self.nStateVars)
        
        self.bftUel.initializeYourself(&self.nodeCoordinates[0])
        
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
        
        with nogil:
            self.bftUel.computeYourself(&U[0], &dU[0],
                                                &Pe[0], 
                                                &Ke[0],
                                                &time[0],
                                                dTime,  
                                                pNewdT[0])
        
    def computeDistributedLoad(self,
                               str loadType,
                               double[::1] P,
                               double[::1] K,
                               int faceID,
                               const double[::1] load,
                               const double[::1] U, 
                               const double[::1] time,
                               double dTime):
        
        self.bftUel.computeDistributedLoad(mapLoadTypes[loadType],
                                    &P[0], 
                                    &K[0], 
                                    faceID,
                                    &load[0],
                                    &U[0], 
                                    &time[0],
                                    dTime)
        
    def computeBodyForce(self, 
            double[::1] P,
            double[::1] K,
           const double[::1] load,
           const double[::1] U, 
           const double[::1] time,
           double dTime):
        
        self.bftUel.computeBodyForce(
                                    &P[0], 
                                    &K[0], 
                                    &load[0],
                                    &U[0], 
                                    &time[0],
                                    dTime)
    def acceptLastState(self,):
        self.stateVars[:] = self.stateVarsTemp
        
    def resetToLastValidState(self,):
        pass
    
    def getResultArray(self, result, gaussPt, getPersistentView=True):    
        """ get the array of a result, possibly as a persistent view which is continiously
        updated by the element """
        cdef string result_ =  result.encode('UTF-8')
        return np.array(  self.getPermanentResultPointer(result_, gaussPt), copy= not getPersistentView)
        
            
    cdef double[::1] getPermanentResultPointer(self, string result, int gaussPt, ):
        """ direct access the the stateVars of the element / underlying material"""
        cdef int resultLength = 0
        cdef double* ptr = self.bftUel.getPermanentResultPointer(result, gaussPt, resultLength)
        return <double[:resultLength]> ( ptr )
    
    def __dealloc__(self):
        del self.bftUel