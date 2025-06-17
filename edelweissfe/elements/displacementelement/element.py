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
from edelweissfe.elements.displacementelement._elementcomputationmatrices import (
    computeBOperator,
    computeJacobian,
    computeNOperator,
)
from edelweissfe.elements.library import elLibrary
from edelweissfe.materials.base.basehyperelasticmaterial import BaseHyperElasticMaterial
from edelweissfe.points.node import Node
from edelweissfe.utils.caseinsensitivedict import CaseInsensitiveDict


class DisplacementElement(BaseElement):
    """This element can be used for EdelweissFE.
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
        if eval(properties["elClass"]) is not DisplacementElement:
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
        self._dStrain = np.zeros([self._nInt, 6])

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

        # initialize the matrices
        self.J = computeJacobian(
            self._xi, self._eta, self._zeta, self._nodesCoordinates, self._nInt, self.nNodes, self.nSpatialDimensions
        )
        self.B = computeBOperator(
            self._xi, self._eta, self._zeta, self._nodesCoordinates, self._nInt, self.nNodes, self.nSpatialDimensions
        )

    def setMaterial(self, material: type):
        """Assign a material.

        Parameters
        ----------
        material
            An initalized instance of a material.
        """

        self.material = material
        if self.planeStrain:  # use 3D
            self._matrixVoigtIndices = np.array([0, 1, 3])
            self._matrixSize = 6
        else:
            self._matrixVoigtIndices = np.arange(self._matrixSize)
        stateVarsSize = 12 + self.material.getNumberOfRequiredStateVars()
        self._dStressdStrain = np.zeros([self._nInt, self._matrixSize, self._matrixSize])
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
            self._isHyperelastic = True
        else:
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

        # assume it's plain strain if it's not given by user
        K = np.reshape(K_, (self._nDof, self._nDof))
        # copy all elements
        self._stateVarsTemp = [self._stateVarsRef[i].copy() for i in range(self._nInt)].copy()
        # strain increment
        self._dStrain[:, self._activeVoigtIndices] = np.array([self.B[i] @ dU for i in range(self._nInt)])
        for i in range(self._nInt):
            # get stress and strain
            stress = self._stateVarsTemp[i][0:6]
            self.material.assignCurrentStateVars(self._stateVarsTemp[i][12:])
            if not self._isHyperelastic:
                # use 3D for 2D planeStrain
                if not self.planeStrain and self.nSpatialDimensions == 2:
                    self.material.computePlaneStress(stress, self._dStressdStrain[i], self._dStrain[i], time, dTime)
                else:
                    self.material.computeStress(stress, self._dStressdStrain[i], self._dStrain[i], time, dTime)
            elif self._isHyperelastic:
                raise Exception("Please use the nonlinear element (displacementtlelement) for hyperelastic materials.")
            # C material tangent
            C = self._dStressdStrain[i][self._matrixVoigtIndices][:, self._matrixVoigtIndices]
            # B operator
            B = self.B[i]
            # Jacobi determinant
            detJ = lin.det(self.J[i])
            # get stiffness matrix for element j in point i
            K += B.T @ C @ B * detJ * self._t * self._weight[i]
            # calculate P
            P -= B.T @ stress[self._activeVoigtIndices] * detJ * self._weight[i] * self._t
            # update strain in stateVars
            self._stateVarsTemp[i][6:12] += self._dStrain[i]

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
