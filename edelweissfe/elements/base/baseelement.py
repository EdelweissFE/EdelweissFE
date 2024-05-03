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
# Created on Fri Aug 5 08:35:06 2022

# @author: Matthias Neuner

"""
Implementing your own finite elements can be done easily by subclassing from 
the abstract base class :class:`~BaseElement`.
"""

from abc import ABC, abstractmethod
import numpy as np
from edelweissfe.points.node import Node
from edelweissfe.timesteppers.timestep import TimeStep


class BaseElement(ABC):
    @property
    @abstractmethod
    def elNumber(self) -> int:
        """The unique number of this element"""

        pass

    @property
    @abstractmethod
    def nNodes(self) -> int:
        """The number of nodes this element requires"""

        pass

    @property
    @abstractmethod
    def nodes(self) -> int:
        """The list of nodes this element holds"""

        pass

    @property
    @abstractmethod
    def nDof(self) -> int:
        """The total number of degrees of freedom this element has"""

        pass

    @property
    @abstractmethod
    def fields(self) -> list[list[str]]:
        """The list of fields per nodes."""

        pass

    @property
    @abstractmethod
    def dofIndicesPermutation(self) -> np.ndarray:
        """The permutation pattern for the residual vector and the stiffness matrix to
        aggregate all entries in order to resemble the defined fields nodewise."""

        pass

    @property
    @abstractmethod
    def ensightType(self) -> str:
        """The shape of the element in Ensight Gold notation."""

        pass

    @property
    @abstractmethod
    def hasMaterial(self) -> str:
        """Flag to check if a material was assigned to this element."""

        pass

    @abstractmethod
    def __init__(self, elementType: str, elNumber: int):
        """Finite elements in EdelweissFE should be derived from this
        base class in order to follow the general interface.

        EdelweissFE expects the layout of the internal and external load vectors, P, PExt, (and the stiffness)
        to be of the form

        .. code-block:: console

            [ node 1 - dofs field 1,
              node 1 - dofs field 2,
              node 1 - ... ,
              node 1 - dofs field n,
              node 2 - dofs field 1,
              ... ,
              node N - dofs field n].

        Parameters
        ----------
        elementType
            A string identifying the requested element formulation.
        elNumber
            A unique integer label used for all kinds of purposes.
        """

        pass

    @abstractmethod
    def setNodes(self, nodes: list[Node]):
        """Assign the nodes to the element.

        Parameters
        ----------
        nodes
            A list of nodes.
        """

        pass

    @abstractmethod
    def setProperties(self, elementProperties: np.ndarray):
        """Assign a set of properties to the element.

        Parameters
        ----------
        elementProperties
            A numpy array containing the element properties.
        """

        pass

    @abstractmethod
    def initializeElement(
        self,
    ):
        """Initalize the element to be ready for computing."""

        pass

    @abstractmethod
    def setMaterial(self, materialName: str, materialProperties: np.ndarray):
        """Assign a material and material properties.

        Parameters
        ----------
        materialName
            The name of the requested material.
        materialProperties
            The numpy array containing the material properties.
        """

    @abstractmethod
    def setInitialCondition(self, stateType: str, values: np.ndarray):
        """Assign initial conditions.

        Parameters
        ----------
        stateType
            The type of initial state.
        values
            The numpy array describing the initial state.
        """

        pass

    @abstractmethod
    def computeDistributedLoad(
        self,
        loadType: str,
        P: np.ndarray,
        K: np.ndarray,
        faceID: int,
        load: np.ndarray,
        U: np.ndarray,
        time: np.ndarray,
        dT: float,
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

        pass

    @abstractmethod
    def computeYourself(
        self,
        P: np.ndarray,
        K: np.ndarray,
        U: np.ndarray,
        dU: np.ndarray,
        time: np.ndarray,
        dT: float,
    ):
        """Evaluate the residual and stiffness for given time, field, and field increment due to a surface load.

        Parameters
        ----------
        P
            The external load vector to be defined.
        K
            The stiffness matrix to be defined.
        U
            The current solution vector.
        dU
            The current solution vector increment.
        time
            Array of step time and total time.
        dTime
            The time increment.
        """

        pass

    @abstractmethod
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

        pass

    @abstractmethod
    def acceptLastState(
        self,
    ):
        """Accept the computed state (in nonlinear iteration schemes)."""

        pass

    @abstractmethod
    def resetToLastValidState(
        self,
    ):
        """Rest to the last valid state."""

        pass

    @abstractmethod
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

        pass

    @abstractmethod
    def getCoordinatesAtCenter(self) -> np.ndarray:
        """Compute the underlying MarmotElement centroid coordinates.

        Returns
        -------
        np.ndarray
            The element's central coordinates.
        """

        pass

    @abstractmethod
    def getCoordinatesAtQuadraturePoints(self) -> np.ndarray:
        """Compute the underlying MarmotElement qp coordinates.

        Returns
        -------
        np.ndarray
            The element's qp coordinates.
        """

        pass

    @abstractmethod
    def getNumberOfQuadraturePoints(self) -> int:
        """Compute the underlying MarmotElement qp coordinates.

        Returns
        -------
        np.ndarray
            The element's qp coordinates.
        """

        pass
