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
#  Daniel Reitmair daniel.reitmair@uibk.ac.at
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

import numpy as np
import numpy.linalg as lin

from edelweissfe.elements.base.baseelement import BaseElement
from edelweissfe.elements.displacementtlelement._elementcomputationmatrices import (
    computeBOperator,
    computeDeformationGradient,
    computeJacobian,
    computeNablaN,
    computeNOperator,
    makeH2D,
)
from edelweissfe.elements.library import elLibrary
from edelweissfe.materials.base.basehyperelasticmaterial import BaseHyperElasticMaterial
from edelweissfe.points.node import Node
from edelweissfe.utils.caseinsensitivedict import CaseInsensitiveDict
from edelweissfe.utils.voigtnotation import doVoigtStrain, undoVoigtStress


def Hgeo(nablaN, S, dim):
    """Computes the unintegrated geometric stiffness matrix for one quadrature point.

    Parameters
    ----------
    nablaN
        The derivative of the shape functions w.r.t the actual coordinates.
    S
        The second Piola-Kirchhoff stress tensor in Voigt notation.
    dim
        Dimension of the domain.

    Returns
    -------
    np.ndarray
        The unintegrated geometric stiffness matrix."""

    S = undoVoigtStress(dim, S)
    Ie = np.eye(dim)
    Hsub = nablaN.T @ S @ nablaN  # [nNodes x nNodes]
    nDof = dim * len(nablaN[0])
    H = np.zeros([nDof, nDof])  # [nDof x nDof]
    for i in range(len(Hsub)):
        for j in range(len(Hsub[0])):
            H[dim * i : dim * (i + 1), dim * j : dim * (j + 1)] = Ie * Hsub[i, j]
    return H


def makeDeformationGradient3D(F, dim):
    """Make the deformation gradient 3D [3 x 3].

    Parameters
    ----------
    F
        The deformation gradient.
    dim
        Dimension of the domain.

    Returns
    -------
    np.ndarray
        The deformation gradient in [3 x 3]."""

    return np.array([[F[0, 0], F[0, 1], 0], [F[1, 0], F[1, 1], 0], [0, 0, 1]]) if dim == 2 else F


class DisplacementTLElement(BaseElement):
    """This total lagrangian element can be used for EdelweissFE.
    The element currently only allows calculations with node forces and given displacements.

    Parameters
    ----------
    elementType
        A string identifying the requested element formulation as shown below.
    elNumber
        A unique integer label used for all kinds of purposes.

    The following types of elements and attributes are currently possible (elementType):

    Elements
    --------
        CPE4
            quadrilateral element with 4 nodes and plane strain.
        CPE8
            quadrilateral element with 8 nodes and plane strain.
        CPS4
            quadrilateral element with 4 nodes and plane stress.
        CPS8
            quadrilateral element with 8 nodes and plane stress.
        C3D8
            hexahedron element with 8 nodes.
        C3D20
            hexahedron element with 20 nodes.

    optional Parameters
    -------------------
    The following attributes are also included in the elType definition:

        R
            reduced integration for element, at the end of elType.
        E
            extended integration for element, at the end of elType.
        N
            (optional) regular integration for element, at the end of elType.

    If R or E is not given by the user, we assume regular increment."""

    @property
    def elNumber(self) -> int:
        """The unique number of this element"""

        return self._elNumber  # return number

    @property
    def nNodes(self) -> int:
        """The number of nodes this element requires"""

        return self._nNodes

    @property
    def nodes(self) -> int:
        """The list of nodes this element holds"""

        return self._nodes

    @property
    def nDof(self) -> int:
        """The total number of degrees of freedom this element has"""

        return self._nDof

    @property
    def fields(self) -> list[list[str]]:
        """The list of fields per nodes."""

        return self._fields

    @property
    def dofIndicesPermutation(self) -> np.ndarray:
        """The permutation pattern for the residual vector and the stiffness matrix to
        aggregate all entries in order to resemble the defined fields nodewise.
        In this case it stays the same because we use the nodes exactly like they are."""

        return self._dofIndices

    @property
    def ensightType(self) -> str:
        """The shape of the element in Ensight Gold notation."""

        return self._ensightType

    @property
    def visualizationNodes(self) -> str:
        """The nodes for visualization."""

        return self._nodes

    @property
    def hasMaterial(self) -> str:
        """Flag to check if a material was assigned to this element."""

        return self._hasMaterial

    def __init__(self, elementType: str, elNumber: int):
        properties = elLibrary[elementType]
        if eval(properties["elClass"]) is not DisplacementTLElement:
            raise Exception("Something went wrong with the element initialization!")
        self._elNumber = elNumber
        self._nNodes = properties["nNodes"]
        self._nDof = properties["nDof"]
        self._dofIndices = properties["dofIndices"]
        self._ensightType = properties["ensightType"]
        self.nSpatialDimensions = properties["nSpatialDimensions"]
        self._nInt = properties["nInt"]
        self._xi = properties["xi"]
        self._eta = properties["eta"]
        self._zeta = properties["zeta"]
        self._weight = properties["w"]
        self._matrixSize = properties["matSize"]
        self._activeVoigtIndices = properties["index"]
        self.planeStrain = properties["plStrain"]
        if self.nSpatialDimensions == 3:
            self._t = 1  # "thickness" for 3D elements
        self._fields = [["displacement"] for i in range(self._nNodes)]
        self._strain = np.zeros([self._nInt, 6])
        self._dStrain = np.zeros([self._nInt, 6])
        # for total Lagrangian formulation
        self._Eold = np.zeros([self._nInt, self.nSpatialDimensions, self.nSpatialDimensions])
        self._E = np.zeros([self._nInt, self.nSpatialDimensions, self.nSpatialDimensions])
        self._F = np.stack(self._nInt * [np.eye(self.nSpatialDimensions)])

    def setNodes(self, nodes: list[Node]):
        """Assign the nodes to the element.

        Parameters
        ----------
        nodes
            A list of nodes.
        """

        self._nodes = nodes
        _nodesCoordinates = np.array([n.coordinates for n in nodes])  # get node coordinates
        self._nodesCoordinates = _nodesCoordinates.transpose()  # nodes given column-wise: x-coordinate - y-coordinate

    def setProperties(self, elementProperties: np.ndarray):
        """Assign a set of properties to the element.

        Parameters
        ----------
        elementProperties
            A numpy array containing the element properties.

        Attributes
        ----------
        thickness
            Thickness of 2D elements.
        """

        if self.nSpatialDimensions == 2:
            self._t = elementProperties[0]  # thickness

    def initializeElement(
        self,
    ):
        """Initalize the element to be ready for computing."""

        self.J = computeJacobian(
            self._xi, self._eta, self._zeta, self._nodesCoordinates, self._nInt, self._nNodes, self.nSpatialDimensions
        )
        self.nablaN = computeNablaN(
            self._xi, self._eta, self._zeta, self.J, self._nInt, self._nNodes, self.nSpatialDimensions
        )

    def setMaterial(self, material: type):
        """Assign a material.

        Parameters
        ----------
        material
            An initalized instance of a material.
        """

        self.material = material
        self._materialProperties = self.material.materialProperties
        if self.planeStrain:  # use 3D
            self._matrixVoigtIndices = np.array([0, 1, 3])
            self._matrixSize = 6
        else:
            self._matrixVoigtIndices = np.arange(self._matrixSize)
        stateVarsSize = 12 + self.material.getNumberOfRequiredStateVars()
        self._hasMaterial = True
        self._stateVarsRef = np.zeros([self._nInt, stateVarsSize])
        self._stateVars = [
            CaseInsensitiveDict(
                {
                    "stress": self._stateVarsRef[i][0:6],
                    "strain": self._stateVarsRef[i][6:12],
                    "materialstate": self._stateVarsRef[i][12:],
                }
            )
            for i in range(self._nInt)
        ]
        self._stateVarsTemp = np.zeros([self._nInt, stateVarsSize])
        if issubclass(type(self.material), BaseHyperElasticMaterial):  # check if material is hyperelastic
            self._dStress_dDeformationGradient = np.zeros([self._nInt, 3, 3, 3, 3])
            self._isHyperelastic = True
        else:
            self._dStress_dStrain = np.zeros([self._nInt, self._matrixSize, self._matrixSize])
            self._isHyperelastic = False

    def setInitialCondition(self, stateType: str, values: np.ndarray):
        """Assign initial conditions.

        Parameters
        ----------
        stateType
            The type of initial state.
        values
            The numpy array describing the initial state.
        """

        raise Exception("Setting an initial condition is not possible with this element provider.")

    def computeDistributedLoad(
        self,
        loadType: str,
        P: np.ndarray,
        K: np.ndarray,
        faceID: int,
        load: np.ndarray,
        U: np.ndarray,
        time: np.ndarray,
        dTime: float,
    ):
        """Evaluate residual and stiffness for given time, field, and field increment due to a surface load.

        Parameters
        ----------
        loadType
            The type of load.
        P
            The external load vector to be defined.
        K
            The stiffness matrix to be defined.
        faceID
            The number of the elements face this load acts on.
        load
            The magnitude (or vector) describing the load.
        U
            The current solution vector.
        time
            Array of step time and total time.
        dTime
            The time increment.
        """

        raise Exception("Applying a distributed load is currently not possible with this element provider.")

    def computeYourself(
        self,
        K_: np.ndarray,
        P: np.ndarray,
        U: np.ndarray,
        dU: np.ndarray,
        time: np.ndarray,
        dTime: float,
    ):
        """Evaluate the residual and stiffness matrix for given time, field, and field increment due to a displacement or load.

        Parameters
        ----------
        P
            The external load vector gets calculated.
        K
            The stiffness matrix gets calculated.
        U
            The current solution vector.
        dU
            The current solution vector increment.
        time
            Array of step time and total time.
        dTime
            The time increment.
        """

        dim = self.nSpatialDimensions
        # assume it's plain strain if it's not given by user
        K = np.reshape(K_, (self._nDof, self._nDof))
        # get current state Vars
        self._stateVarsTemp = [self._stateVarsRef[i].copy() for i in range(self._nInt)].copy()
        # compute the deformation gradient
        self._F = computeDeformationGradient(U, self.nablaN, self._nInt, self._nNodes, dim)
        if not self._isHyperelastic:
            B = computeBOperator(self._F, self.nablaN, self._nInt, self._nNodes, dim)
        for i in range(self._nInt):
            detJ = lin.det(self.J[i])
            # get stress (PK2) and strain
            stress = self._stateVarsTemp[i][0:6]
            self.material.assignCurrentStateVars(self._stateVarsTemp[i][12:])
            H = self._F[i] - np.eye(dim)
            self._E[i] = 1 / 2 * (H + H.T + H.T @ H)
            invF = lin.inv(self._F[i])
            F = makeDeformationGradient3D(self._F[i], dim)
            if self._isHyperelastic:
                NAi = np.zeros([self._nNodes, 3])
                _nablaN = np.zeros([3, self._nNodes])
                NAi[:, :dim] = self.nablaN[i].T @ invF
                _nablaN[:dim] = self.nablaN[i]
                self._strain[i, self._activeVoigtIndices] = doVoigtStrain(dim, self._E[i])
                # use 3D for 2D planeStrain
                if not self.planeStrain and dim == 2:
                    self.material.computePlaneKirchhoff(stress, self._dStress_dDeformationGradient[i], F, time, dTime)
                    T = undoVoigtStress(2, stress)
                else:
                    invF = makeDeformationGradient3D(invF, dim)
                    self.material.computeKirchhoff(stress, self._dStress_dDeformationGradient[i], F, time, dTime)
                    T = undoVoigtStress(3, stress)
                PK1 = invF @ T
                # update strain in stateVars
                self._stateVarsTemp[i][6:12] = self._strain[i]
                Hk = makeH2D(
                    np.einsum("ai,ijkl,bl->ajbk", NAi, self._dStress_dDeformationGradient[i], _nablaN.T)
                    - np.einsum("ak,bi,ij->ajbk", NAi, NAi, T),
                    dim,
                )
                # compute inner forces
                P -= (self.nablaN[i].T @ PK1[:dim, :dim]).flatten() * detJ * self._t * self._weight[i]
            else:  # for non-hyperelastic materials
                self._dStrain[i, self._activeVoigtIndices] = doVoigtStrain(dim, self._E[i] - self._Eold[i])
                if not self.planeStrain and dim == 2:
                    self.material.computePlaneStress(stress, self._dStress_dStrain[i], self._dStrain[i], time, dTime)
                else:
                    self.material.computeStress(stress, self._dStress_dStrain[i], self._dStrain[i], time, dTime)
                self._stateVarsTemp[i][6:12] += self._dStrain[i]
                Cm = self._dStress_dStrain[i][self._matrixVoigtIndices][:, self._matrixVoigtIndices]
                Hk = B[i].T @ Cm @ B[i] + Hgeo(self.nablaN[i], stress, dim)
                # compute inner forces
                P -= B[i].T @ stress[self._matrixVoigtIndices] * detJ * self._weight[i] * self._t
            # calculate complete stiffness matrix
            K += Hk * detJ * self._t * self._weight[i]

    def computeBodyForce(
        self, P: np.ndarray, K: np.ndarray, load: np.ndarray, U: np.ndarray, time: np.ndarray, dTime: float
    ):
        """Evaluate residual and stiffness for given time, field, and field increment due to a body force load.

        Parameters
        ----------
        P
            The external load vector to be defined.
        K
            The stiffness matrix to be defined.
        load
            The magnitude (or vector) describing the load.
        U
            The current solution vector.
        time
            Array of step time and total time.
        dTime
            The time increment.
        """

        N = computeNOperator(self._xi, self._eta, self._zeta, self._nInt, self.nNodes, self.nSpatialDimensions)
        for i in range(self._nInt):
            P += np.outer(N[i], load).flatten() * lin.det(self.J[i]) * self._t * self._weight[i]

    def acceptLastState(
        self,
    ):
        """Accept the computed state (in nonlinear iteration schemes)."""

        # copy every array in array (complete copying)
        self._stateVarsRef[:] = [self._stateVarsTemp[i][:] for i in range(self._nInt)]
        self._Eold = self._E.copy()

    def resetToLastValidState(
        self,
    ):
        """Reset to the last valid state."""

    def getResultArray(self, result: str, quadraturePoint: int, getPersistentView: bool = True) -> np.ndarray:
        """Get the array of a result, possibly as a persistent view which is continiously
        updated by the element.

        Parameters
        ----------
        result
            The name of the result.
        quadraturePoint
            The number of the quadrature point.
        getPersistentView
            If true, the returned array should be continiously updated by the element.

        Returns
        -------
        np.ndarray
            The result.
        """

        try:
            return self._stateVars[quadraturePoint][result]
        except KeyError:  # result in material
            self.material.assignCurrentStateVars(self._stateVarsRef[quadraturePoint][12:])
            return self.material.getResult(result)

    def getCoordinatesAtCenter(self) -> np.ndarray:
        """Compute the underlying MarmotElement centroid coordinates.

        Returns
        -------
        np.ndarray
            The element's central coordinates.
        """

        x = self._nodesCoordinates
        return np.average(x, axis=1)

    def getNumberOfQuadraturePoints(self) -> int:
        """Get the number of Quadrature points the element has.

        Returns
        -------
        nInt
            The number of Quadrature points.
        """

        return self._nInt

    def getCoordinatesAtQuadraturePoints(self) -> np.ndarray:
        """Compute the underlying MarmotElement qp coordinates.

        Returns
        -------
        np.ndarray
            The element's qp coordinates.
        """

        N = computeNOperator(self._xi, self._eta, self._zeta, self._nInt, self.nNodes, self.nSpatialDimensions)
        return self._nodesCoordinates @ N
