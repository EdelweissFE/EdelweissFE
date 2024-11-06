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
# Created on Thu Apr 27 08:35:06 2017

# @author: matthias

import numpy as np

cimport cython
cimport libcpp.cast
cimport numpy as np

cimport edelweissfe.elements.marmotelement.element

from edelweissfe.utils.exceptions import CutbackRequest

from libc.stdlib cimport free, malloc
from libcpp.memory cimport allocator, make_unique, unique_ptr

from edelweissfe.elements.base.baseelement import BaseElement

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
    # cdef classes cannot subclass. Hence we do not subclass from the BaseElement,
    # but still we follow the interface for compatiblity.

    def __init__(self, elementType, elNumber):
        """This element serves as a wrapper for MarmotElements.

        For the documentation of MarmotElements, please refer to `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_.

        Parameters
        ----------
        elementType
            The Marmot element which should be represented, e.g., CPE4.
        elNumber
            The number of the element."""

        self._elNumber = elNumber
        self._elType = elementType

        self._nNodes                         = self.marmotElement.getNNodes()
        self._nSpatialDimensions             = self.marmotElement.getNSpatialDimensions()
        self._nDof                           = self.marmotElement.getNDofPerElement()

        cdef vector[vector[string]] fields  = self.marmotElement.getNodeFields()
        self._fields                         = [ [ s.decode('utf-8')  for s in n  ] for n in fields ]

        cdef vector[int] permutationPattern = self.marmotElement.getDofIndicesPermutationPattern()
        self._dofIndicesPermutation          = np.asarray(permutationPattern)

        self._ensightType                    = self.marmotElement.getElementShape().decode('utf-8')

        self._hasMaterial = False

    def __cinit__(self, elementType, elNumber):
        """This C-level method is responsible for actually creating the MarmotElement.

        Parameters
        ----------
        elementType
            The Marmot element which should be represented, e.g., CPE4.
        elNumber
            The number of the element."""

        try:
            self.marmotElement = MarmotElementFactory.createElement(MarmotElementFactory.getElementCodeFromName( elementType.upper().encode('utf-8')), self._elNumber)
        except IndexError:
            raise NotImplementedError("Marmot element {:} not found in library.".format(elementType))

    @property
    def elNumber(self):
        return self._elNumber

    @property
    def nSpatialDimensions(self):
       return self._nSpatialDimensions

    @property
    def elType(self):
        return self._elType

    @property
    def nodes(self):
        return self._nodes

    @property
    def nNodes(self):
        return self._nNodes

    @property
    def nDof(self):
        return self._nDof

    @property
    def fields(self):
        return self._fields

    @property
    def dofIndicesPermutation(self):
        return self._dofIndicesPermutation

    @property
    def ensightType(self):
        return self._ensightType

    @property
    def visualizationNodes(self):
        return self._nodes

    @property
    def hasMaterial(self):
        return self._hasMaterial

    def setNodes (self, nodes):
        """Assign the nodes coordinates to the underyling MarmotElement"""

        self._nodes = nodes
        self.nodeCoordinates = np.concatenate([ node.coordinates for node in nodes])
        self.marmotElement.assignNodeCoordinates(&self.nodeCoordinates[0])

    def setProperties(self, elementProperties):
        """Assign a set of properties to the underyling MarmotElement"""

        self._elementProperties = elementProperties

        self.marmotElement.assignProperty(
                ElementProperties(
                        &self._elementProperties[0],
                        self._elementProperties.shape[0] ) )

    def initializeElement(self, ):
        """Let the underlying MarmotElement initialize itself"""
        self.marmotElement.initializeYourself()

    def setMaterial(self, materialName, materialProperties):
        """Assign a material and material properties to the underlying MarmotElement.
        Furthermore, create two sets of state vars:

            - the actual set,
            - and a temporary set for backup in nonlinear iteration schemes.
        """
        self._materialProperties =  materialProperties
        try:
            self.marmotElement.assignProperty(
                    MarmotMaterialSection(
                            MarmotMaterialFactory.getMaterialCodeFromName(
                                    materialName.upper().encode('UTF-8')),
                            &self._materialProperties[0],
                            self._materialProperties.shape[0] ) )
        except IndexError:
            raise NotImplementedError("Marmot material {:} not found in library.".format(materialName))

        self.nStateVars =           self.marmotElement.getNumberOfRequiredStateVars()

        self._stateVars =            np.zeros(self.nStateVars)
        self._stateVarsTemp =        np.zeros(self.nStateVars)

        self.marmotElement.assignStateVars(&self._stateVarsTemp[0], self.nStateVars)

        self._hasMaterial = True

    cpdef void _initializeStateVarsTemp(self, ) nogil:
        self._stateVarsTemp[:] = self._stateVars

    def setInitialCondition(self,
                            stateType,
                            const double[::1] values):
        """Assign initial conditions to the underlying Marmot element"""

        if not self._hasMaterial:
            raise Exception("Element {:} has no material assigned!".format(self._elNumber))

        self._initializeStateVarsTemp()
        self.marmotElement.setInitialConditions(mapStateTypes[stateType], &values[0])
        self.acceptLastState()

    cpdef void computeYourself(self,
                         double[::1] Ke,
                         double[::1] Pe,
                         const double[::1] U,
                         const double[::1] dU,
                         const double[::1] time,
                         double dTime, ) nogil except *:
        """Evaluate residual and stiffness for given time, field, and field increment."""

        if not self._hasMaterial:
            raise Exception("Element {:} has no material assigned!".format(self._elNumber))

        cdef double pNewDT
        with nogil:
            self._initializeStateVarsTemp()

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
        """Evaluate residual and stiffness for given time, field, and field increment due to a surface load."""

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
        """Evaluate residual and stiffness for given time, field, and field increment due to a volume load."""

        self.marmotElement.computeBodyForce(
                                    &P[0],
                                    &K[0],
                                    &load[0],
                                    &U[0],
                                    &time[0],
                                    dTime)
    def acceptLastState(self,):
        """Accept the computed state (in nonlinear iteration schemes)."""

        self._stateVars[:] = self._stateVarsTemp

    def resetToLastValidState(self,):
        """Reset to the last valid state."""


    def getResultArray(self, result, quadraturePoint, getPersistentView=True):
        """Get the array of a result, possibly as a persistent view which is continiously
        updated by the underlying MarmotElement."""

        if not self._hasMaterial:
            raise Exception("Element {:} has no material assigned!".format(self._elNumber))

        cdef string result_ =  result.encode('UTF-8')
        return np.array(  self.getStateView(result_, quadraturePoint), copy= not getPersistentView)

    cdef double[::1] getStateView(self, string result, int quadraturePoint, ):
        """Directly access the state vars of the underlying MarmotElement"""

        if not self._hasMaterial:
            raise Exception("Element {:} has no material assigned!".format(self._elNumber))

        cdef StateView res = self.marmotElement.getStateView(result, quadraturePoint )

        return <double[:res.stateSize]> ( res.stateLocation )

    def getCoordinatesAtCenter(self):
        """Compute the underlying MarmotElement centroid coordinates."""

        return np.asarray ( self.marmotElement.getCoordinatesAtCenter() )

    def getCoordinatesAtQuadraturePoints(self):
        """Compute the underlying MarmotElement qp coordinates."""

        return np.asarray ( self.marmotElement.getCoordinatesAtQuadraturePoints() )

    def getNumberOfQuadraturePoints(self):
        """Compute the underlying MarmotElement qp coordinates."""

        return self.marmotElement.getNumberOfQuadraturePoints()

    def __dealloc__(self):
        del self.marmotElement
