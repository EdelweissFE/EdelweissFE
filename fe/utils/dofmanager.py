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
# Created on Tue Dec 18 09:18:25 2018

# @author: Matthias Neuner
"""
This module contains important classes for describing the global equation system by means of a sparse system.
"""

from collections import OrderedDict
from fe.config.phenomena import getFieldSize, phenomena
from fe.variables.node import Node
import numpy as np


class VIJSystemMatrix(np.ndarray):
    """
    This class represents the V Vector of VIJ triple (sparse matrix in COO format),
    which

      * also contains the I and J vectors as class members,
      * allows to directly access (contiguous read and write) access of each entity via the [] operator

    Parameters
    ----------
    nDof
        The size of the system.
    I
        The I vector for the VIJ triple.
    J
        The J vector for the VIJ triple.
    entitiesInVIJ
        A dictionary containing the indices of an entitiy in the value vector.
    """

    def __new__(cls, nDof: int, I: np.ndarray, J: np.ndarray, entitiesInVIJ: dict):

        obj = np.zeros_like(I, dtype=float).view(cls)

        obj.nDof = nDof
        obj.I = I
        obj.J = J
        obj.entitiesInVIJ = entitiesInVIJ

        return obj

    def __getitem__(self, key):
        try:
            idxInVIJ = self.entitiesInVIJ[key]
            return super().__getitem__(slice(idxInVIJ, idxInVIJ + key.nDof**2))
        except:
            return super().__getitem__(key)


class DofVector(np.ndarray):
    """
    This class represents a Dof Vector, which also has knowledge of each entities (elements, constraints) location within.
    The [] operator allows to access (non-contigouos read, write) at each entities location

    Parameters
    ----------
    nDof
        The size of the system.
    entitiesInVIJ
        A dictionary containing the indices of an entitiy in the value vector.
    """

    def __new__(cls, nDof: int, entitiesInDofVector: dict):
        obj = np.zeros(nDof, dtype=float).view(cls)
        obj.entitiesInDofVector = entitiesInDofVector

        return obj

    def __getitem__(self, key):
        try:
            return super().__getitem__(self.entitiesInDofVector[key])
        except:
            return super().__getitem__(key)

    def __setitem__(self, key, value):
        try:
            return super().__setitem__(self.entitiesInDofVector[key], value)
        except:
            return super().__setitem__(key, value)


class DofManager:
    """
    The DofManager

     * analyzes the domain (nodes and constraints),
     * collects information about the necessary structure of the degrees of freedom
     * handles the active fields on each node
     * counts the accumulated number of associated elements on each dof (for the Abaqus like convergence test)
     * supplies the framework with DofVectors and VIJSystemMatrices

    Parameters
    ----------
    model
        A dictionary containing the model tree.
    """

    def __init__(self, model: dict):

        self.model = model

        self.activateFieldsOnNodes()

        (
            self.nDof,
            self.indicesOfFieldsInDofVector,
            self.indicesOfScalarVariablesInDofVector,
        ) = self.initializeDofVectorStructure()

        (self.nAccumulatedNodalFluxes, self.nAccumulatedNodalFluxesFieldwise) = self.countNodalFluxes()

        (
            self.accumulatedElementNDof,
            self.accumulatedConstraintNDof,
            self.sizeVIJ,
            self.largestNumberOfElNDof,
        ) = self.countAccumulatedEntityDofs()

        self.entitiesInDofVector = self.locateEntitiesInDofVector()

        (self.I, self.J, self.entitiesInVIJ) = self.initializeVIJPattern()

    def activateFieldsOnNodes(
        self,
    ):
        """Activate all fields on nodes, which are required in the analysis."""

        model = self.model

        for element in model["elements"].values():
            for node, nodeFields in zip(element.nodes, element.fields):
                for field in nodeFields:
                    node.fields[field] = True

        for constraint in model["constraints"].values():
            for node, nodeFields in zip(constraint.nodes, constraint.fieldsOnNodes):
                for field in nodeFields:
                    node.fields[field] = True

    def initializeDofVectorStructure(self) -> tuple[int, dict[str, np.ndarray]]:
        """Loop over all nodes to generate the global field-dof indices.

        Returns
        -------
        tuple
            output is a tuple of:
             * number of total DOFS
             * dictionary of fields and indices:
                * field
                * indices
        """

        nodes = self.model["nodes"]
        domainSize = self.model["domainSize"]
        # constraints = self.model["constraints"]

        indicesOfFieldsInDofVector = OrderedDict()
        currentIndexInDofVector = 0

        # blockwise assembly
        for field in phenomena.keys():
            fieldSize = getFieldSize(field, domainSize)
            indexList = []

            for node in nodes.values():
                if field in node.fields:
                    node.fields[field] = [i + currentIndexInDofVector for i in range(fieldSize)]
                    indexList += node.fields[field]
                    currentIndexInDofVector += fieldSize

            # if we have dofs of this field at all
            if indexList:
                indicesOfFieldsInDofVector[field] = indexList

        self.indexToNodeMapping = {
            index: node for node in nodes.values() for field in node.fields.values() for index in field
        }

        indicesOfScalarVariablesInDofVector = []
        for scalarVariable in self.model["scalarVariables"]:
            indicesOfScalarVariablesInDofVector.append(currentIndexInDofVector)
            scalarVariable.index = currentIndexInDofVector
            currentIndexInDofVector += 1

        nDof = currentIndexInDofVector

        return nDof, indicesOfFieldsInDofVector, indicesOfScalarVariablesInDofVector

    def countNodalFluxes(self) -> tuple[int, dict[str, int]]:
        """For the VIJ (COO) system matrix and the Abaqus like convergence test,
        the number of dofs 'entity-wise' is needed:
        = Σ_(elements+constraints) Σ_nodes ( nDof (field) ).

        Returns
        -------
        tuple
            - Number of accumulated fluxes in total
            - Number of accumulated fluxes per field:
                - Field
                - Number of accumulated fluxes
        """

        indicesOfFieldsInDofVector = self.indicesOfFieldsInDofVector

        elements = self.model["elements"]
        constraints = self.model["constraints"]

        accumulatedNumberOfFieldFluxes = {}
        accumulatedNodalFluxesTotal = 0

        for field in indicesOfFieldsInDofVector.keys():
            accumulatedNumberOfFieldFluxes[field] = np.sum(
                [len(node.fields[field]) for el in elements.values() for node in (el.nodes) if field in node.fields]
                + [
                    len(node.fields[field])
                    for constraint in constraints.values()
                    for node in (constraint.nodes)
                    if field in node.fields
                ]
            )

            accumulatedNodalFluxesTotal += accumulatedNumberOfFieldFluxes[field]

        return (
            accumulatedNodalFluxesTotal,
            accumulatedNumberOfFieldFluxes,
        )

    def countAccumulatedEntityDofs(self) -> tuple[int, int, int, int]:
        """Generates some auxiliary information,
        which may be required by some modules of EdelweissFE.

        Returns
        -------
        tuple[int,int,int,int]
            The tuple of
                - number of elemental degrees of freedom,
                - number of constraint degrees of freedom,
                - size of system matrix,
                - largest number of dofs on a element
        """

        elements = self.model["elements"]
        constraints = self.model["constraints"]

        sizeVIJ = 0
        accumulatedElementNDof = 0
        largestNumberOfElNDof = 0

        for el in elements.values():
            accumulatedElementNDof += el.nDof
            sizeVIJ += el.nDof**2

            largestNumberOfElNDof = max(el.nDof, largestNumberOfElNDof)

        nNDofAccumulatedConstraints = 0
        for constraint in constraints.values():
            nNDofAccumulatedConstraints += constraint.nDof

            sizeVIJ += constraint.nDof**2

        return accumulatedElementNDof, nNDofAccumulatedConstraints, sizeVIJ, largestNumberOfElNDof

    def locateEntitiesInDofVector(
        self,
    ) -> dict:
        """Creates a dictionary containing the location (indices) of each entity (elements, constraints)
        within the DofVector structure.

        Returns
        -------
        dict
            A dictionary containing the location mapping.
        """

        elements = self.model["elements"]
        constraints = self.model["constraints"]
        entitiesInDofVector = {}

        for el in elements.values():
            destList = np.asarray(
                [
                    i
                    for iNode, node in enumerate(el.nodes)  # for each node of the element..
                    for nodeField in el.fields[iNode]  # for each field of this node
                    for i in node.fields[nodeField]
                ]
            )  # the index in the global system

            entitiesInDofVector[el] = destList[el.dofIndicesPermutation]

        for constraint in constraints.values():
            # destList = constraint.globalDofIndices
            destList = np.asarray(
                [
                    i
                    for iNode, node in enumerate(constraint.nodes)  # for each node of the constraint
                    for nodeField in constraint.fieldsOnNodes[iNode]  # for each field of this node
                    for i in node.fields[nodeField]
                ]
                + [v.index for v in constraint.scalarVariables]
            )
            entitiesInDofVector[constraint] = destList

        return entitiesInDofVector

    def initializeVIJPattern(
        self,
    ) -> tuple[np.ndarray, np.ndarray, dict]:
        """Generate the IJ pattern for VIJ (COO) system matrices.

        Returns
        -------
        tuple
             - I vector
             - J vector
             - the entities to system matrix entry mapping.
        """

        elements = self.model["elements"]
        constraints = self.model["constraints"]

        entitiesInVIJ = {}
        entitiesInDofVector = self.entitiesInDofVector

        sizeVIJ = self.sizeVIJ

        I = np.zeros(sizeVIJ, dtype=int)
        J = np.zeros(sizeVIJ, dtype=int)
        idxInVIJ = 0

        for el in elements.values():
            destList = entitiesInDofVector[el]

            entitiesInVIJ[el] = idxInVIJ

            # looks like black magic, but it's an efficient way to generate all indices of Ke in K:
            elDofLocations = np.tile(destList, (destList.shape[0], 1))
            I[idxInVIJ : idxInVIJ + el.nDof**2] = elDofLocations.ravel()
            J[idxInVIJ : idxInVIJ + el.nDof**2] = elDofLocations.ravel("F")
            idxInVIJ += el.nDof**2

        for constraint in constraints.values():
            destList = entitiesInDofVector[constraint]

            entitiesInVIJ[constraint] = idxInVIJ

            constraintDofLocations = np.tile(destList, (destList.shape[0], 1))
            I[idxInVIJ : idxInVIJ + constraint.nDof**2] = constraintDofLocations.ravel()
            J[idxInVIJ : idxInVIJ + constraint.nDof**2] = constraintDofLocations.ravel("F")
            idxInVIJ += constraint.nDof**2

        return I, J, entitiesInVIJ

    def constructVIJSystemMatrix(
        self,
    ) -> VIJSystemMatrix:
        """Construct a VIJ (COO) Sparse System matrix object, which also has knowledge about
        the location of each entity.

        Returns
        -------
        VIJSystemMatrix
            The system Matrix.
        """

        nDof = self.nDof
        I = self.I
        J = self.J
        entitiesInVIJ = self.entitiesInVIJ

        return VIJSystemMatrix(nDof, I, J, entitiesInVIJ)

    def constructDofVector(
        self,
    ) -> DofVector:
        """Construct a vector with size=nDof and which has knowledge about
        the location of each entity.

        Returns
        -------
        DofVector
            A DofVector.
        """

        nDof = self.nDof
        entitiesInDofVector = self.entitiesInDofVector

        return DofVector(nDof, entitiesInDofVector)

    def getNodeForIndexInDofVector(self, index: int) -> Node:
        """Find the node for a given index in the equuation system."""

        return self.indexToNodeMapping[index]
