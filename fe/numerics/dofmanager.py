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
from fe.points.node import Node
from fe.fields.nodefield import NodeField
from fe.models.femodel import FEModel
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

    def __init__(self, model: FEModel):
        self.initializeFromModel(model)

    def _initializeDofVectorStructure(self, model: FEModel) -> tuple[int, dict[str, np.ndarray]]:
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

        indicesOfFieldsInDofVector = dict()
        indicesOfNodeFieldVariablesInDofVector = dict()
        currentIndexInDofVector = 0

        for fieldName, field in model.nodeFields.items():
            nextIndexInDofVector = currentIndexInDofVector + field.dimension * len(field.nodes)
            indicesOfFieldsInDofVector[fieldName] = slice(currentIndexInDofVector, nextIndexInDofVector)

            indicesOfNodeFieldVariablesInDofVector |= {
                n.fields[fieldName]: np.arange(
                    currentIndexInDofVector + i * field.dimension,
                    currentIndexInDofVector + i * field.dimension + field.dimension,
                    dtype=int,
                )
                for i, n in enumerate(field.nodes)
            }
            currentIndexInDofVector = nextIndexInDofVector

        indicesOfScalarVariablesInDofVector = {
            scalarVariable: currentIndexInDofVector + i
            for i, scalarVariable in enumerate(model.scalarVariables.values())
        }

        currentIndexInDofVector += len(indicesOfScalarVariablesInDofVector)
        nDof = currentIndexInDofVector

        return (
            nDof,
            indicesOfNodeFieldVariablesInDofVector,
            indicesOfFieldsInDofVector,
            indicesOfScalarVariablesInDofVector,
        )

    def _determineIndexToNodeMapping(self, model):
        indexToNodeMapping = {
            index: node
            for node in model.nodes.values()
            for f in node.fields.values()
            for index in self.idcsOfFieldVariablesInDofVector[f]
            if f in self.idcsOfFieldVariablesInDofVector
        }

        return indexToNodeMapping

    def _countNodalFluxes(self, model) -> tuple[int, dict[str, int]]:
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

        indicesOfFieldsInDofVector = self.idcsOfFieldsInDofVector

        elements = model.elements
        constraints = model.constraints

        accumulatedNumberOfFieldFluxes = {}
        accumulatedNodalFluxesTotal = 0

        for field in indicesOfFieldsInDofVector.keys():
            accumulatedNumberOfFieldFluxes[field] = np.sum(
                [
                    len(self.idcsOfFieldVariablesInDofVector[node.fields[field]])
                    for el in elements.values()
                    for node in (el.nodes)
                    if field in node.fields
                ]
                + [
                    len(self.idcsOfFieldVariablesInDofVector[node.fields[field]])
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

    def _countAccumulatedEntityDofs(self, model) -> tuple[int, int, int, int]:
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

        elements = model.elements
        constraints = model.constraints

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

    def _locateElementsInDofVector(self, model) -> dict:
        """Creates a dictionary containing the location (indices) of each entity (elements, constraints)
        within the DofVector structure.

        Returns
        -------
        dict
            A dictionary containing the location mapping.
        """

        elements = model.elements
        idcsOfElementsInDofVector = {}

        for el in elements.values():
            destList = np.hstack(
                [
                    self.idcsOfFieldVariablesInDofVector[node.fields[nodeField]]
                    for iNode, node in enumerate(el.nodes)  # for each node of the element..
                    for nodeField in el.fields[iNode]  # for each field of this node
                ]
            )  # the index in the global system

            idcsOfElementsInDofVector[el] = destList[el.dofIndicesPermutation]

        return idcsOfElementsInDofVector

    def _locateConstraintsInDofVector(self, model) -> dict:
        """Creates a dictionary containing the location (indices) of each entity (elements, constraints)
        within the DofVector structure.

        Returns
        -------
        dict
            A dictionary containing the location mapping.
        """

        constraints = model.constraints
        idcsOfConstraintsInDofVector = {}

        for constraint in constraints.values():
            destList = np.hstack(
                [
                    self.idcsOfFieldVariablesInDofVector[node.fields[nodeField]]
                    for iNode, node in enumerate(constraint.nodes)  # for each node of the constraint
                    for nodeField in constraint.fieldsOnNodes[iNode]  # for each field of this node
                ]
                + [self.idcsOfScalarVariablesInDofVector[v] for v in constraint.scalarVariables]
            )
            idcsOfConstraintsInDofVector[constraint] = destList

        return idcsOfConstraintsInDofVector

    def _locateFieldsOnNodeSetFieldsInDofVector(self, model) -> dict:
        """Creates a dictionary containing the location (indices) of each entity (elements, constraints)
        within the DofVector structure.

        Returns
        -------
        dict
            A dictionary containing the location mapping.
        """

        nodeSets = model.nodeSets
        nodeSetFieldsInDofVector = {}

        for field in self.idcsOfFieldsInDofVector:
            nodeSetFieldsInDofVector[field] = dict()

            for nSetName, nSet in nodeSets.items():
                nodeSetFieldsInDofVector[field][nSetName] = np.array(
                    [self.idcsOfFieldVariablesInDofVector[node.fields[field]] for node in nSet if field in node.fields]
                ).flatten()

        return nodeSetFieldsInDofVector

    def _initializeVIJPattern(
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

        entitiesInVIJ = {}
        entitiesInDofVector = self.idcsOfHigherOrderEntitiesInDofVector

        sizeVIJ = self.sizeVIJ

        I = np.zeros(sizeVIJ, dtype=int)
        J = np.zeros(sizeVIJ, dtype=int)
        idxInVIJ = 0

        for entity, entityIdcsInDofVector in self.idcsOfHigherOrderEntitiesInDofVector.items():
            entitiesInVIJ[entity] = idxInVIJ

            nDofEntity = len(entityIdcsInDofVector)

            # looks like black magic, but it's an efficient way to generate all indices of Ke in K:
            VIJLocations = np.tile(entityIdcsInDofVector, (nDofEntity, 1))
            I[idxInVIJ : idxInVIJ + nDofEntity**2] = VIJLocations.flatten()
            J[idxInVIJ : idxInVIJ + nDofEntity**2] = VIJLocations.flatten("F")
            idxInVIJ += nDofEntity**2

        return I, J, entitiesInVIJ

    def getIndicesOfAllVariablesInDofVector(self):
        """
        Get list of indices of all lower order entitties (variables).

        Returns
        -------
        list
            The indices.
        """
        return self.idcsOfFieldVariablesInDofVector | self.idcsOfScalarVariablesInDofVector

    def getIndicesOfAllHigherOrderEntitiesInDofVector(self):
        """
        Get list of indices of all higher order entitties.

        Returns
        -------
        list
            The indices.
        """

        return self.idcsOfElementsInDofVector | self.idcsOfConstraintsInDofVector

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

        return DofVector(self.nDof, self.idcsInDofVector)

    def getNodeForIndexInDofVector(self, index: int) -> Node:
        """Find the node for a given index in the equuation system.

        Parameters
        ----------
        index
            The index in the DofVector.

        Returns
        -------
        Node
            The attached Node.
        """

        return self.indexToNodeMapping[index]

    def writeDofVectorToNodeField(self, dofVector, nodeField, resultName):
        """Write the current values of an entire NodeField from the respective locations in a given DofVector.


        Parameters
        ----------
        dofVector
            The source DofVector.
        nodeField
            The NodeField to get the updated values.
        resultname
            The name of the value entries held by the NodeField.

        Returns
        -------
        NodeField
            The updated NodeField.
        """

        if not resultName in nodeField.values:
            nodeField.createFieldValueEntry(resultName)

        nodeField.values[resultName][:] = dofVector[self.idcsOfFieldsInDofVector[nodeField.name]].reshape(
            (-1, nodeField.dimension)
        )

        return nodeField

    def writeNodeFieldToDofVector(self, dofVector: DofVector, nodeField: NodeField, resultName: str):
        """Write the current values of an entire NodeField to the respective locations in a given DofVector.


        Parameters
        ----------
        dofVector
            The result DofVector.
        nodeField
            The NodeField holding the values.
        resultname
            The name of the value entries held by the NodeField.

        Returns
        -------
        DofVector
            The DofVector.
        """

        dofVector[self.idcsOfFieldsInDofVector[nodeField.name]] = nodeField.values[resultName].flatten()

        return dofVector

    def initializeFromModel(self, model: FEModel):
        """Initialize the manged DofVector and VIJ-format system matrix structures from a given FEModel tree.


        Parameters
        ----------
        model
            The model tree.
        """

        (
            self.nDof,
            self.idcsOfFieldVariablesInDofVector,
            self.idcsOfFieldsInDofVector,
            self.idcsOfScalarVariablesInDofVector,
        ) = self._initializeDofVectorStructure(model)

        self.indexToNodeMapping = self._determineIndexToNodeMapping(model)

        (self.nAccumulatedNodalFluxes, self.nAccumulatedNodalFluxesFieldwise) = self._countNodalFluxes(model)

        (
            self.accumulatedElementNDof,
            self.accumulatedConstraintNDof,
            self.sizeVIJ,
            self.largestNumberOfElNDof,
        ) = self._countAccumulatedEntityDofs(model)

        self.idcsOfFieldsOnNodeSetsInDofVector = self._locateFieldsOnNodeSetFieldsInDofVector(model)
        self.idcsOfElementsInDofVector = self._locateElementsInDofVector(model)
        self.idcsOfConstraintsInDofVector = self._locateConstraintsInDofVector(model)

        self.idcsOfVariablesInDofVector = self.getIndicesOfAllVariablesInDofVector()
        self.idcsOfHigherOrderEntitiesInDofVector = self.getIndicesOfAllHigherOrderEntitiesInDofVector()

        self.idcsInDofVector = self.idcsOfVariablesInDofVector | self.idcsOfHigherOrderEntitiesInDofVector

        (self.I, self.J, self.entitiesInVIJ) = self._initializeVIJPattern()
