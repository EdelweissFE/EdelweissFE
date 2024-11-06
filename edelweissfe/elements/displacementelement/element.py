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

import numpy as np
import numpy.linalg as lin

from edelweissfe.elements.base.baseelement import BaseElement
from edelweissfe.elements.displacementelement.elementmatrices import (
    computeBOperator,
    computeJacobian,
    computeNOperator,
)
from edelweissfe.points.node import Node
from edelweissfe.utils.caseinsensitivedict import CaseInsensitiveDict

# element parameters 2D
w1 = (5 / 9) ** 2
w2 = (5 / 9) * (8 / 9)
w3 = (8 / 9) ** 2
# element parameters 3D
s8 = 1 / np.sqrt(3) * np.array([-1, 1, 1, -1])  # get s
t8 = 1 / np.sqrt(3) * np.array([-1, -1, 1, 1])  # get z
s20 = np.array([-1, 0, 1])  # get s
t20 = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1])  # get t
w1H = (5 / 9) ** 3
w2H = (5 / 9) ** 2 * (8 / 9)
w3H = (5 / 9) * (8 / 9) ** 2
w4H = (8 / 9) ** 3
wI = np.array([w1H, w2H, w1H])
wII = np.array([w2H, w3H, w2H])
# Element library
elLibrary = CaseInsensitiveDict(
    Quad4=dict(
        nNodes=4,
        nDof=8,
        dofIndices=np.arange(0, 8),
        ensightType="quad4",
        nSpatialDimensions=2,
        nInt=4,
        xi=1 / np.sqrt(3) * np.array([1, 1, -1, -1]),
        eta=1 / np.sqrt(3) * np.array([1, -1, -1, 1]),
        zeta=None,
        w=np.ones(4),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=True,
    ),
    Quad4R=dict(
        nNodes=4,
        nDof=8,
        dofIndices=np.arange(0, 8),
        ensightType="quad4",
        nSpatialDimensions=2,
        nInt=1,
        xi=np.array([0]),
        eta=np.array([0]),
        zeta=None,
        w=np.array([4]),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=True,
    ),
    Quad4E=dict(
        nNodes=4,
        nDof=8,
        dofIndices=np.arange(0, 8),
        ensightType="quad4",
        nSpatialDimensions=2,
        nInt=9,
        xi=np.sqrt(0.6) * np.array([0, -1, -1, 1, 1, -1, 0, 1, 0]),
        eta=np.sqrt(0.6) * np.array([0, -1, 1, 1, -1, 0, 1, 0, -1]),
        zeta=None,
        w=np.array([w3, w1, w1, w1, w1, w2, w2, w2, w2]),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=True,
    ),
    Quad8=dict(
        nNodes=8,
        nDof=16,
        dofIndices=np.arange(0, 16),
        ensightType="quad8",
        nSpatialDimensions=2,
        nInt=9,
        xi=np.sqrt(0.6) * np.array([0, -1, -1, 1, 1, -1, 0, 1, 0]),
        eta=np.sqrt(0.6) * np.array([0, -1, 1, 1, -1, 0, 1, 0, -1]),
        zeta=None,
        w=np.array([w3, w1, w1, w1, w1, w2, w2, w2, w2]),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=True,
    ),
    Quad8R=dict(
        nNodes=8,
        nDof=16,
        dofIndices=np.arange(0, 16),
        ensightType="quad8",
        nSpatialDimensions=2,
        nInt=4,
        xi=1 / np.sqrt(3) * np.array([1, 1, -1, -1]),
        eta=1 / np.sqrt(3) * np.array([1, -1, -1, 1]),
        zeta=None,
        w=np.ones(4),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=True,
    ),
    Quad4PS=dict(
        nNodes=4,
        nDof=8,
        dofIndices=np.arange(0, 8),
        ensightType="quad4",
        nSpatialDimensions=2,
        nInt=4,
        xi=1 / np.sqrt(3) * np.array([1, 1, -1, -1]),
        eta=1 / np.sqrt(3) * np.array([1, -1, -1, 1]),
        zeta=None,
        w=np.ones(4),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=False,
    ),
    Quad4RPS=dict(
        nNodes=4,
        nDof=8,
        dofIndices=np.arange(0, 8),
        ensightType="quad4",
        nSpatialDimensions=2,
        nInt=1,
        xi=np.array([0]),
        eta=np.array([0]),
        zeta=None,
        w=np.array([4]),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=False,
    ),
    Quad4EPS=dict(
        nNodes=4,
        nDof=8,
        dofIndices=np.arange(0, 8),
        ensightType="quad4",
        nSpatialDimensions=2,
        nInt=9,
        xi=np.sqrt(0.6) * np.array([0, -1, -1, 1, 1, -1, 0, 1, 0]),
        eta=np.sqrt(0.6) * np.array([0, -1, 1, 1, -1, 0, 1, 0, -1]),
        zeta=None,
        w=np.array([w3, w1, w1, w1, w1, w2, w2, w2, w2]),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=False,
    ),
    Quad8PS=dict(
        nNodes=8,
        nDof=16,
        dofIndices=np.arange(0, 16),
        ensightType="quad8",
        nSpatialDimensions=2,
        nInt=9,
        xi=np.sqrt(0.6) * np.array([0, -1, -1, 1, 1, -1, 0, 1, 0]),
        eta=np.sqrt(0.6) * np.array([0, -1, 1, 1, -1, 0, 1, 0, -1]),
        zeta=None,
        w=np.array([w3, w1, w1, w1, w1, w2, w2, w2, w2]),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=False,
    ),
    Quad8RPS=dict(
        nNodes=8,
        nDof=16,
        dofIndices=np.arange(0, 16),
        ensightType="quad8",
        nSpatialDimensions=2,
        nInt=4,
        xi=1 / np.sqrt(3) * np.array([1, 1, -1, -1]),
        eta=1 / np.sqrt(3) * np.array([1, -1, -1, 1]),
        zeta=None,
        w=np.ones(4),
        matSize=3,
        index=np.array([0, 1, 3]),
        plStrain=False,
    ),
    Hexa8=dict(
        nNodes=8,
        nDof=24,
        dofIndices=np.arange(0, 24),
        ensightType="hexa8",
        nSpatialDimensions=3,
        nInt=8,
        xi=1 / np.sqrt(3) * np.hstack([-np.ones(4), np.ones(4)]),
        eta=np.hstack([t8, t8]),
        zeta=np.hstack([s8, s8]),
        w=np.ones(8),
        matSize=6,
        index=np.arange(6),
        plStrain=None,
    ),
    Hexa8R=dict(
        nNodes=8,
        nDof=24,
        dofIndices=np.arange(0, 24),
        ensightType="hexa8",
        nSpatialDimensions=3,
        nInt=1,
        xi=np.array([0]),
        eta=np.array([0]),
        zeta=np.array([0]),
        w=np.array([8]),
        matSize=6,
        index=np.arange(6),
        plStrain=None,
    ),
    Hexa8E=dict(
        nNodes=8,
        nDof=24,
        dofIndices=np.arange(0, 24),
        ensightType="hexa8",
        nSpatialDimensions=3,
        nInt=27,
        xi=np.sqrt(0.6) * np.hstack([-np.ones(9), np.zeros(9), np.ones(9)]),
        eta=np.sqrt(0.6) * np.hstack([t20, t20, t20]),
        zeta=np.sqrt(0.6) * np.hstack([s20 for i in range(9)]),
        w=np.hstack([wI, wII, wI, wII, np.array([w3H, w4H, w3H]), wII, wI, wII, wI]),
        matSize=6,
        index=np.arange(6),
        plStrain=None,
    ),
    Hexa20=dict(
        nNodes=20,
        nDof=60,
        dofIndices=np.arange(0, 60),
        ensightType="hexa20",
        nSpatialDimensions=3,
        nInt=27,
        xi=np.sqrt(0.6) * np.hstack([-np.ones(9), np.zeros(9), np.ones(9)]),
        eta=np.sqrt(0.6) * np.hstack([t20, t20, t20]),
        zeta=np.sqrt(0.6) * np.hstack([s20 for i in range(9)]),
        w=np.hstack([wI, wII, wI, wII, np.array([w3H, w4H, w3H]), wII, wI, wII, wI]),
        matSize=6,
        index=np.arange(6),
        plStrain=None,
    ),
    Hexa20R=dict(
        nNodes=20,
        nDof=60,
        dofIndices=np.arange(0, 60),
        ensightType="hexa20",
        nSpatialDimensions=3,
        nInt=8,
        xi=1 / np.sqrt(3) * np.hstack([-np.ones(4), np.ones(4)]),
        eta=np.hstack([t8, t8]),
        zeta=np.hstack([s8, s8]),
        w=np.ones(8),
        matSize=6,
        index=np.arange(6),
        plStrain=None,
    ),
)
# add variations
elLibrary.update(
    {
        "Quad4PE": elLibrary["Quad4"],
        "Quad4RPE": elLibrary["Quad4R"],
        "Quad4EPE": elLibrary["Quad4E"],
        "Quad8PE": elLibrary["Quad8"],
        "Quad8RPE": elLibrary["Quad8R"],
    }
)


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
        Quad4
            quadrilateral element with 4 nodes.
        Quad8
            quadrilateral element with 8 nodes.
        Hexa8
            hexahedron element with 8 nodes.

    optional Parameters
    -------------------
    The following attributes are also included in the elementtype definition:

        R
            reduced integration for element, in elementtype[5].
        E
            extended integration for element, in elementtype[5].
        PE
            use plane strain for 2D elements, in elementtype[6:8] or [5:7].
        PS
            use plane stress for 2D elements, in elementtype[6:8] or [5:7].

    If R or E is not given by the user, we assume regular increment.
    If PE or PS is not given by the user, we assume PE."""

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
        self.elementtype = elementType[0].upper() + elementType[1:5].lower() + elementType[5:].upper()
        self._elNumber = elNumber
        try:
            if len(self.elementtype) > 5 and self.elementtype[5].lower() == "n":  # normal integration
                self.elementtype = self.elementtype.replace("N", "").replace("n", "")
            properties = elLibrary[self.elementtype]
        except KeyError:
            raise Exception("This element type doesn't exist.")
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

        if self.elementtype[0] == "Q":
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
            # use 3D for 2D planeStrain
            if not self.planeStrain and self.nSpatialDimensions == 2:
                self.material.computePlaneStress(stress, self._dStressdStrain[i], self._dStrain[i], time, dTime)
            else:
                self.material.computeStress(stress, self._dStressdStrain[i], self._dStrain[i], time, dTime)
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

        return self._stateVars[quadraturePoint][result]

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
