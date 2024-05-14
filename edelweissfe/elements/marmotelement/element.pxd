#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  ---------------------------------------------------------------------
#
#  _____    _      _              _         _____ _____
# | ____|__| | ___| |_      _____(_)___ ___|  ___| ____|
# |  _| / _` |/ _ \ \ \ /\ / / _ \ / __/ __| |_  |  _|
# | |__| (_| |  __/ |\ V  V /  __/ \__ \__ \  _| | |___
# |_____\__,_|\___|_| \_/\_/ \___|_|___/___/_|   |_____|
#
#
#  Unit of Strength of Materials and Structural Analysis
#  University of Innsbruck,
#  2017 - today
#
#  Matthias Neuner matthias.neuner@uibk.ac.at
#
#  This file is part of EdelweissFE.
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  The full text of the license can be found in the file LICENSE.md at
#  the top level directory of EdelweissFE.
#  ---------------------------------------------------------------------
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""
cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np


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
    cdef cppclass MarmotMaterialFactory:
        @staticmethod
        int getMaterialCodeFromName(const string& materialName) except +IndexError

    cdef cppclass MarmotElementFactory:
        @staticmethod
        int getElementCodeFromName(const string& elementName) except +IndexError
        @staticmethod
        MarmotElement* createElement(int elementCode, int noEl,) except +ValueError

cdef extern from "Marmot/MarmotElementProperty.h":
    cdef cppclass MarmotElementProperty nogil:
        pass

    cdef cppclass MarmotMaterialSection(MarmotElementProperty) nogil:
        MarmotMaterialSection(int materialCode, const double* _materialProperties, int nMaterialProperties)

    cdef cppclass ElementProperties(MarmotElementProperty) nogil:
        ElementProperties(const double* _elementProperties, int nElementProperties)

cdef extern from "Marmot/MarmotUtils.h":
    cdef struct StateView:
        double *stateLocation
        int stateSize


cdef extern from "Marmot/MarmotElement.h":
    cdef cppclass MarmotElement nogil:

        int getNumberOfRequiredStateVars()

        void assignStateVars(double *_stateVars, int nStateVars)

        void assignProperty( const MarmotElementProperty& property )

        void assignProperty( const MarmotMaterialSection& property ) except +ValueError

        void assignNodeCoordinates(const double* elementCoordinates)

        void initializeYourself()

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

        StateView getStateView(const string& stateName, int gaussPt)

        vector[vector[string]] getNodeFields()

        vector[int] getDofIndicesPermutationPattern()

        string getElementShape()

        int getNNodes()

        int getNSpatialDimensions()

        int getNDofPerElement()

        vector[double] getCoordinatesAtCenter()

        vector[vector[double]] getCoordinatesAtQuadraturePoints()

        int getNumberOfQuadraturePoints()

cdef class MarmotElementWrapper:

    cdef MarmotElement* marmotElement
    cdef list _nodes,
    cdef int _elNumber,
    cdef str _elType,
    cdef int _nNodes, _nDof, _nSpatialDimensions
    cdef list _fields
    cdef str _ensightType
    cdef int _hasMaterial
    cdef np.ndarray _dofIndicesPermutation

    cdef public double[::1] _stateVars, nodeCoordinates
    cdef public double[::1] _elementProperties, _stateVarsTemp , _materialProperties
    cdef int nStateVars

    cdef double[::1] getStateView(self, string stateName, int gaussPt)

    # nogil methods are already declared here:

    cpdef void _initializeStateVarsTemp(self, ) nogil

    cpdef void computeYourself(self,
                     double[::1] Ke,
                     double[::1] Pe,
                     const double[::1] U,
                     const double[::1] dU,
                     const double[::1] time,
                     double dTime, ) nogil except *
