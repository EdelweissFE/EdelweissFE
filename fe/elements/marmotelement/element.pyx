#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

import numpy as np
cimport numpy as np
cimport fe.elements.marmotelement.element 
cimport libcpp.cast
cimport cython

from fe.utils.exceptions import CutbackRequest

from libcpp.memory cimport unique_ptr, allocator, make_unique

from libc.stdlib cimport malloc, free

mapLoadTypes={
        'pressure' : DistributedLoadTypes.Pressure,
        'surface torsion' : DistributedLoadTypes.SurfaceTorsion,
        'surface traction' : DistributedLoadTypes.SurfaceTraction
     }

mapStateTypes={
        'geostatic stress' : StateTypes.GeostaticStress,
        'sdvini' : StateTypes.MarmotMaterialStateVars,
        'initialize material': StateTypes.MarmotMaterialInitialization
     }
    
@cython.final # no subclassing -> cpdef with nogil possible
cdef class MarmotElementWrapper:
    
    def __cinit__(self, elementType, nodes, elNumber):
        self.nodes = nodes
        if self.nodes:
            self.nodeCoordinates = np.concatenate([ node.coordinates for node in nodes])
            
        self.elNumber = elNumber
        
        self.marmotElement = MarmotElementFactory.createElement( MarmotElementFactory.getElementCodeFromName( elementType.upper().encode('utf-8')), self.elNumber)
        
        self.nNodes                         = self.marmotElement.getNNodes()
        
        self.nDof                           = self.marmotElement.getNDofPerElement()
        
        cdef vector[vector[string]] fields  = self.marmotElement.getNodeFields()
        self.fields                         = [ [ s.decode('utf-8')  for s in n  ] for n in fields ]
        
        cdef vector[int] permutationPattern = self.marmotElement.getDofIndicesPermutationPattern()
        self.dofIndicesPermutation          = np.asarray(permutationPattern)
        
        self.ensightType                    = self.marmotElement.getElementShape().decode('utf-8')
        
    def setProperties(self, elementProperties, materialName, materialProperties):
        
        self.elementProperties =    elementProperties
        self.materialProperties =   materialProperties
        
        self.marmotElement.assignProperty( 
                ElementProperties(
                        &self.elementProperties[0],
                        self.elementProperties.shape[0] ) )
        
        self.marmotElement.assignProperty(
                MarmotMaterialSection(
                        MarmotMaterialFactory.getMaterialCodeFromName(
                                materialName.upper().encode('UTF-8')), 
                        &self.materialProperties[0],
                        self.materialProperties.shape[0] ) )
        
        self.nStateVars =           self.marmotElement.getNumberOfRequiredStateVars()
        
        self.stateVars =            np.zeros(self.nStateVars)
        self.stateVarsTemp =        np.zeros(self.nStateVars)
        
        self.marmotElement.assignStateVars(&self.stateVarsTemp[0], self.nStateVars)
        
        self.marmotElement.initializeYourself(&self.nodeCoordinates[0])
        
    cpdef void initializeStateVarsTemp(self, ) nogil:
        self.stateVarsTemp[:] = self.stateVars
        
    def setInitialCondition(self, 
                            stateType,
                            const double[::1] values):
        
        self.initializeStateVarsTemp()
        self.marmotElement.setInitialConditions(mapStateTypes[stateType], &values[0])
        self.acceptLastState()

    cpdef void computeYourself(self, 
                         double[::1] Ke, 
                         double[::1] Pe, 
                         const double[::1] U, 
                         const double[::1] dU, 
                         const double[::1] time, 
                         double dTime, ) nogil except *:
        cdef double pNewDT 
        with nogil:
            self.initializeStateVarsTemp()
            
            pNewDT = 1e36
            
            self.marmotElement.computeYourself(&U[0], &dU[0],
                                                &Pe[0], 
                                                &Ke[0],
                                                &time[0],
                                                dTime,  
                                                pNewDT)
            if pNewDT < 1.0:
                raise CutbackRequest("Element {:} requests for a cutback!".format(self.elNumber), pNewDT)
        
    def computeDistributedLoad(self,
                               str loadType,
                               double[::1] P,
                               double[::1] K,
                               int faceID,
                               const double[::1] load,
                               const double[::1] U, 
                               const double[::1] time,
                               double dTime):
        
        self.marmotElement.computeDistributedLoad(mapLoadTypes[loadType],
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
        
        self.marmotElement.computeBodyForce(
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
        cdef PermanentResultLocation res = self.marmotElement.getPermanentResultPointer(result, gaussPt )

        # TODO: check if we do not 'need' to cast away the constness:
        return <double[:res.resultLength]> ( res.resultLocation)
    
    def __dealloc__(self):
        del self.marmotElement
