#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Marmot/MarmotElement.h" namespace "MarmotElement":
    cdef enum StateTypes:
        Sigma11,
        Sigma22,
        Sigma33,
        HydrostaticStress,
        GeostaticStress,
        MarmotMaterialStateVars,
        MarmotMaterialInitialization
        
    cdef enum DistributedLoadTypes:
        Pressure
        SurfaceTraction
        SurfaceTorsion
        
cdef extern from "Marmot/Marmot.h" namespace "MarmotLibrary" nogil:
    enum MaterialCode: pass
    enum ElementCode: pass

    cdef cppclass MarmotMaterialFactory:
        @staticmethod
        MaterialCode getMaterialCodeFromName(const string& materialName) except +IndexError
    
    cdef cppclass MarmotElementFactory:
        @staticmethod
        ElementCode  getElementCodeFromName(const string& elementName) except +IndexError
        @staticmethod
        MarmotElement* createElement(ElementCode elementCode, int noEl,) except +ValueError
                       
cdef extern from "Marmot/MarmotElementProperty.h":
    cdef cppclass MarmotElementProperty nogil:
        pass
    
    cdef cppclass MarmotMaterialSection(MarmotElementProperty) nogil:
        MarmotMaterialSection(int materialCode, const double* materialProperties, int nMaterialProperties)
        
    cdef cppclass ElementProperties(MarmotElementProperty) nogil:
        ElementProperties(const double* elementProperties, int nElementProperties)

cdef extern from "Marmot/MarmotUtils.h":
    cdef struct PermanentResultLocation:
        const double *resultLocation
        int resultLength


cdef extern from "Marmot/MarmotElement.h":
    cdef cppclass MarmotElement nogil:
        
        int getNumberOfRequiredStateVars()

        void assignStateVars(double *stateVars, int nStateVars)
        
        void assignProperty( const MarmotElementProperty& property ) 

        void assignProperty( const MarmotMaterialSection& property ) except +ValueError

        void initializeYourself(const double* elementCoordinates)

        void computeYourself( const double* QTotal,
                                            const double* dQ,
                                            double* Pe,
                                            double* Ke,
                                            const double* time,
                                            double dT,
                                            double& pNewdT,) except +ValueError
        
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
        
        PermanentResultLocation getPermanentResultPointer(const string& resultName, int gaussPt)
        
        vector[vector[string]] getNodeFields()

        vector[int] getDofIndicesPermutationPattern()

        string getElementShape()
        
        int getNNodes()
    
        int getNDofPerElement()
        
cdef class MarmotElementWrapper:
    
    cdef MarmotElement* marmotElement
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
