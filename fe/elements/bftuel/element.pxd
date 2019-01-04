#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""
from libcpp.string cimport string
from libcpp.vector cimport vector

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
        

cdef extern from "userLibrary.h" namespace "userLibrary" nogil:
    enum MaterialCode: pass
    enum ElementCode: pass

    MaterialCode getMaterialCodeFromName(const string& materialName) except +ValueError
    ElementCode  getElementCodeFromName(const string& elementName) except +ValueError
    
    BftUel* UelFactory(int elementCode, 
                       int noEl,
                       ) except +ValueError
                       
cdef extern from "bftUelProperty.h":
    cdef cppclass BftUelProperty nogil:
        pass
    
    cdef cppclass BftMaterialSection(BftUelProperty) nogil:
        BftMaterialSection(int materialCode, const double* materialProperties, int nMaterialProperties)
        
    cdef cppclass ElementProperties(BftUelProperty) nogil:
        ElementProperties(const double* elementProperties, int nElementProperties)

cdef extern from "bftUel.h":
    cdef cppclass BftUel nogil:
        
        int getNumberOfRequiredStateVars()

        void assignStateVars(double *stateVars, int nStateVars)
        
        void assignProperty( const BftUelProperty& property ) 

        void assignProperty( const BftMaterialSection& property ) 

        void initializeYourself(const double* elementCoordinates)

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
                                double* K, 
                                int faceID, 
                                const double* load,
                                const double* QTotal,
                                const double* time,
                                double dT)
        
        void computeBodyForce(
                        double* P, 
                        double* K, 
                        const double* load,
                        const double* QTotal,
                        const double* time,
                        double dT)
        
        double* getPermanentResultPointer(const string& resultName, int gaussPt, int& resultLength)
        
        vector[vector[string]] getNodeFields()

        vector[int] getDofIndicesPermutationPattern()

        string getElementShape()
        
        int getNNodes()
    
        int getNDofPerElement()
        
cdef class BftUelWrapper:
    
    cdef BftUel* bftUel
    cdef public nodes, 
    cdef public int elNumber, 
    cdef public int nNodes, nDof
    cdef public list fields
    cdef public str ensightType
    cdef readonly dofIndicesPermutation
    
    cdef public double[::1] stateVars, nodeCoordinates
    cdef double[::1] elementProperties, stateVarsTemp , materialProperties
    cdef int nStateVars
    cdef double[::1] getPermanentResultPointer(self, string result, int gaussPt)
    
    # nogil methods are already declared here:
    
    cpdef void initializeStateVarsTemp(self, ) nogil

    cpdef void computeYourself(self, 
                     double[::1] Ke, 
                     double[::1] Pe, 
                     const double[::1] U, 
                     const double[::1] dU, 
                     const double[::1] time, 
                     double dTime, ) nogil except *
