#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

#from fe.materials.umatlibrary cimport pUmatType
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
        
    cdef enum PropertyTypes:
            ElementProperties,
            BftMaterial

cdef extern from "userLibrary.h" namespace "userLibrary" nogil:
    enum MaterialCode: pass
    enum ElementCode: pass

    MaterialCode getMaterialCodeFromName(const string& materialName) except +ValueError
    ElementCode  getElementCodeFromName(const string& elementName) except +ValueError
    
    BftUel* UelFactory(int elementCode, 
                       int noEl,
                       ) except +ValueError

cdef extern from "bftUel.h":
    cdef cppclass BftUel nogil:
        
        int getNumberOfRequiredStateVars()

        void assignStateVars(double *stateVars, int nStateVars)
        
        void assignProperty(PropertyTypes property, int propertyInfo, const double* propertyValues, int nProperties)

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
                                int faceID, 
                                const double* load,
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
    cdef public int nNodes, nDofPerEl
    cdef public list fields
    cdef public str ensightType
    cdef readonly dofIndicesPermutation
    
    cdef public double[::1] stateVars, nodeCoordinates
    cdef double[::1] elementProperties, stateVarsTemp , materialProperties
    cdef int nStateVars
    cdef double[::1] getPermanentResultPointer(self, string result, int gaussPt)
