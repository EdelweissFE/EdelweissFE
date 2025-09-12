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

import numpy as np

from edelweissfe.config.phenomena import phenomena
from edelweissfe.fields.nodefield import NodeField
from edelweissfe.nodecouplingentity.base.nodecouplingentity import (
    BaseNodeCouplingEntity,
)
from edelweissfe.points.node import Node
from edelweissfe.sets.nodeset import NodeSet
from edelweissfe.variables.scalarvariable import ScalarVariable


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

    def __new__(cls, nDof: int, I: np.ndarray, J: np.ndarray, entitiesInVIJ: dict):  # noqa: E741
        obj = np.zeros_like(I, dtype=float).view(cls)

        obj.nDof = nDof
        obj.I = I  # noqa: E741
        obj.J = J
        obj.entitiesInVIJ = entitiesInVIJ

        return obj

    def __getitem__(self, key):
        try:
            idxInVIJ = self.entitiesInVIJ[key]
            return super().__getitem__(slice(idxInVIJ, idxInVIJ + key.nDof**2))
        except Exception:
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
        except Exception:
            return super().__getitem__(key)

    def __setitem__(self, key, value):
        try:
            return super().__setitem__(self.entitiesInDofVector[key], value)
        except Exception:
            return super().__setitem__(key, value)

    def copy(self):
        """Create a copy of the DofVector, including the entitiesInDofVector mapping.

        Returns
        -------
        DofVector
            The copied DofVector.
        """
        newDofVector = super().copy().view(DofVector)
        newDofVector.entitiesInDofVector = self.entitiesInDofVector.copy()
        return newDofVector


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
    nodeFields
        The list of NodeFields which should be represented in the DofVector structure.
    scalarVariables
        The list of ScalarVariables which should be represented in the DofVector structure.
    elements
        The list of Elements for which a map to the respective indices should be created.
    constraints
        The list of Constraints for which a map to the respective indices should be created.
    nodeSets
        The list of NodeSets for which a map to the respective indices should be created.
    """

    def __init__(
        self,
        nodeFields: list[NodeField],
        scalarVariables: list[ScalarVariable] = [],
        elements: list[BaseNodeCouplingEntity] = [],
        constraints: list[BaseNodeCouplingEntity] = [],
        nodeSets: list[NodeSet] = [],
        initializeVIJPattern: bool = True,
    ):

        self.nDof = int()  #: The total number of degrees of freedom (and size of the DofVector)
        self.fields = list()  #: The list of fields which can be found in the Dofvector
        self.idcsOfFieldsInDofVector = dict()  #: The dictionary mapping the field names to all indices in the DofVector
        self.idcsOfScalarVariablesInDofVector = (
            dict()
        )  #: The dictionary mapping the scalar variables to all indices in the DofVector

        self.idcsOfFieldVariablesInDofVector = (
            dict()
        )  #: The dictionary mapping the nodal field variables to the indices in the DofVector
        self.idcsOfNodeFieldsInDofVector = (
            dict()
        )  #: The dictionary mapping the a complete NodeField to the all its indices in the DofVector
        self.indexToHostObjectMapping = (
            dict()
        )  #: The reverse dictionary mapping an index to the Host (e.g., a Node) holding the index's FieldvVariable
        self.accumulatedElementNDof = (
            int()
        )  #: The accumulated number of element dofs (= the sum of all element vector sizes)
        self.largestNumberOfElNDof = int()  #: The size of the largest of all element dof vectors
        self.accumulatedConstraintNDof = (
            int()
        )  #: The accumulated number of constraint dofs (= the sum of all constraint vector sizes)
        self.largestNumberOfConstraintNDof = int()  #: The size of the largest of all constraint dof vectors

        self.idcsOfFieldsOnNodeSetsInDofVector = (
            dict()
        )  #: The dictionary mapping for each field a NodeSet to the respective indices in the DofVector
        self.idcsOfElementsInDofVector = dict()  #: The dictionary mapping an element to it's indices in the DofVector
        self.idcsOfConstraintsInDofVector = (
            dict()
        )  #: The dictionary mapping a constraint to it's indices in the DofVector

        self._sizeVIJ = (
            int()
        )  #: The number of nonzero entries of the system matrix, resulting from the dense higher order entities contributions

        # initialization:

        self._determineDofsAndTheirIndices(nodeFields, scalarVariables)
        self._gatherInformationAboutEntities(elements, constraints, nodeSets)

        self.idcsOfFieldsInDofVector = self.idcsOfNodeFieldsInDofVector
        self.fields = self.idcsOfFieldsInDofVector.keys()

        self.nAccumulatedNodalFluxesFieldwise = self._computeAccumulatedNodalFluxesFieldWise(self.fields)

        self.idcsOfBasicVariablesInDofVector = (
            self.idcsOfFieldVariablesInDofVector | self.idcsOfScalarVariablesInDofVector
        )
        self.idcsOfHigherOrderEntitiesInDofVector = self.idcsOfElementsInDofVector | self.idcsOfConstraintsInDofVector

        self._sizeVIJ = self._accumulatedElementVIJSize + self._accumulatedConstraintVIJSize

        if initializeVIJPattern:
            (self.I, self.J, self.idcsOfHigherOrderEntitiesInVIJ) = self._initializeVIJPattern()

    def _determineDofsAndTheirIndices(self, nodeFields: list, scalarVariables: list):

        (
            self._nDofNodeFields,
            self.idcsOfFieldVariablesInDofVector,
            self.idcsOfNodeFieldsInDofVector,
        ) = self._reserveSpaceForNodeFields(self.nDof, nodeFields)

        self.nDof += self._nDofNodeFields

        self.indexToHostObjectMapping |= self._determineIndexToNodeMap()

        (
            self._nDofScalarVariables,
            self.idcsOfScalarVariablesInDofVector,
        ) = self._reserveSpaceForScalarVariables(self.nDof, scalarVariables)

        self.nDof += self._nDofScalarVariables

    def _gatherInformationAboutEntities(self, elements, constraints, nodeSets):

        (
            self.accumulatedElementNDof,
            self._accumulatedElementVIJSize,
            self._nAccumulatedNodalFluxesFieldwiseFromElements,
            self.largestNumberOfElNDof,
        ) = self._gatherElementsInformation(elements)

        (
            self.accumulatedConstraintNDof,
            self._accumulatedConstraintVIJSize,
            self._nAccumulatedNodalFluxesFieldwiseFromConstraints,
            self.largestNumberOfConstraintNDof,
        ) = self._gatherConstraintsInformation(constraints)

        self.idcsOfFieldsOnNodeSetsInDofVector = self._locateFieldsOnNodeSetsInDofVector(nodeSets)
        self.idcsOfElementsInDofVector = self._locateNodeCouplingEntitiesInDofVector(elements)
        self.idcsOfConstraintsInDofVector = self._locateConstraintsInDofVector(constraints)

    def _reserveSpaceForNodeFields(
        self,
        idxStart: int,
        nodeFields: list[NodeField],
    ) -> tuple[int, dict[str, np.ndarray]]:
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

        idcsOfFieldsInDofVector = dict()
        idcsOfNodeFieldVariablesInDofVector = dict()
        currentIdxInDofVector = idxStart

        for nodeField in nodeFields:
            nextIdxInDofVector = currentIdxInDofVector + nodeField.dimension * len(nodeField.nodes)
            idcsOfFieldsInDofVector[nodeField.name] = slice(currentIdxInDofVector, nextIdxInDofVector)

            idcsOfNodeFieldVariablesInDofVector |= {
                n.fields[nodeField.name]: np.arange(
                    currentIdxInDofVector + i * nodeField.dimension,
                    currentIdxInDofVector + i * nodeField.dimension + nodeField.dimension,
                    dtype=int,
                )
                for i, n in enumerate(nodeField.nodes)
            }
            currentIdxInDofVector = nextIdxInDofVector

        nDof = currentIdxInDofVector

        return (
            nDof,
            idcsOfNodeFieldVariablesInDofVector,
            idcsOfFieldsInDofVector,
        )

    def _reserveSpaceForScalarVariables(
        self, idxStart: int, scalarVariables: list
    ) -> tuple[int, dict[str, np.ndarray]]:
        """Loop over all ScalarVariables to generate their global indices.

        Returns
        -------
        tuple
            output is a tuple of:
             * number of total DOFS
             * dictionary of fields and indices:
                * field
                * indices
        """

        currentIdxInDofVector = idxStart

        idcsOfScalarVariablesInDofVector = {
            scalarVariable: currentIdxInDofVector + i for i, scalarVariable in enumerate(scalarVariables)
        }

        nDof = len(idcsOfScalarVariablesInDofVector)

        return (
            nDof,
            idcsOfScalarVariablesInDofVector,
        )

    def _determineIndexToNodeMap(
        self,
    ) -> dict[int, Node]:
        """Determine the map from each index (associated with a FieldVariable)
        in the DofVector to the corresponding attached Node oject.
        Useful for determining, e.g., the Node associated with a residual outlier in nonlinear
        simulations.

        Returns
        -------
        dict[int, Node]
            The dictionary containing the map from index of a FieldVariable (component) in the DofVector to the respective Node instance.
        """
        indexToNodeMapping = dict()

        for (
            fieldVariable,
            fieldVariableIndices,
        ) in self.idcsOfFieldVariablesInDofVector.items():
            for index in fieldVariableIndices:
                indexToNodeMapping[index] = fieldVariable.node

        return indexToNodeMapping

    def _computeAccumulatedNodalFluxesFieldWise(self, fields: list) -> dict:
        """For the VIJ (COO) system matrix and the Abaqus like convergence test,
        the number of dofs 'entity-wise' is needed:
        = Σ_(elements+constraints) Σ_nodes ( nDof (field) ).

        Parameters
        ----------
        fields:
            The list of fields for which the accumulated nodal fluxes should be computed.

        Returns
        -------
        dict
            Number of accumulated fluxes per field:
                - Field
                - Number of accumulated fluxes
        """

        nAccumulatedNodalFluxesFieldwise = {}

        for field in fields:
            nAccumulatedNodalFluxesFieldwise[field] = self._nAccumulatedNodalFluxesFieldwiseFromElements.get(
                field, 0
            ) + self._nAccumulatedNodalFluxesFieldwiseFromConstraints.get(field, 0)

        return nAccumulatedNodalFluxesFieldwise

    def _gatherElementsInformation(self, entities: list) -> tuple[int, int, int, int]:
        """Generates some auxiliary information,
        which may be required by some modules of EdelweissFE.

        Parameters
        ----------
        entities
           The list of entities, for which the information is gathered.

        Returns
        -------
        tuple[int,int]
            The tuple of
                - number of accumulated elemental degrees of freedom.
                - number of accumulated system matrix sizes.
                - the number of  acummulated fluxes Σ_entities Σ_nodes ( nDof (field) ) for Abaqus-like convergence tests.
                - largest occuring number of dofs on any element.
        """
        accumulatedEntityNDof = 0
        accumulatedEntityVIJSize = 0
        largestNumberOfAnyEntitityDof = 0

        nAccumulatedFluxesFieldwise = dict.fromkeys(phenomena.keys(), 0)

        for e in entities:
            accumulatedEntityNDof += e.nDof
            accumulatedEntityVIJSize += e.nDof**2

            for node in e.nodes:
                for field, fv in node.fields.items():
                    nAccumulatedFluxesFieldwise[field] += len(self.idcsOfFieldVariablesInDofVector[fv])

            largestNumberOfAnyEntitityDof = max(e.nDof, largestNumberOfAnyEntitityDof)

        return (
            accumulatedEntityNDof,
            accumulatedEntityVIJSize,
            nAccumulatedFluxesFieldwise,
            largestNumberOfAnyEntitityDof,
        )

    def _gatherConstraintsInformation(self, entities: list) -> tuple[int, int, int, int]:
        """Generates some auxiliary information,
        which may be required by some modules of EdelweissFE.

        Parameters
        ----------
        entities
           The list of entities, for which the information is gathered.

        Returns
        -------
        tuple[int,int]
            The tuple of
                - number of accumulated elemental degrees of freedom.
                - number of accumulated system matrix sizes.
                - the number of  acummulated fluxes Σ_entities Σ_nodes ( nDof (field) ) for Abaqus-like convergence tests.
                - largest occuring number of dofs on any element.
        """

        return self._gatherElementsInformation(entities)

    # def _analyzeVIJPattern(self,):

    def _locateNodeCouplingEntitiesInDofVector(self, entities: list) -> dict:
        """Creates a dictionary containing the location (indices) of each entity (elements, constraints)
        within the DofVector structure.

        Returns
        -------
        dict
            A dictionary containing the location mapping.
        """

        idcsOfElementsInDofVector = {}

        for ent in entities:
            destList = np.hstack(
                [
                    self.idcsOfFieldVariablesInDofVector[node.fields[nodeField]]
                    for iNode, node in enumerate(ent.nodes)  # for each node of the element..
                    for nodeField in ent.fields[iNode]  # for each field of this node
                ]
            )  # the index in the global system

            if ent.dofIndicesPermutation is not None:
                idcsOfElementsInDofVector[ent] = destList[ent.dofIndicesPermutation]
            else:
                idcsOfElementsInDofVector[ent] = destList

        return idcsOfElementsInDofVector

    def _locateConstraintsInDofVector(self, constraints: list) -> dict:
        """Creates a dictionary containing the location (indices) of each entity (elements, constraints)
        within the DofVector structure.

        Returns
        -------
        dict
            A dictionary containing the location mapping.
        """

        constraints = constraints
        idcsOfConstraintsInDofVector = {}

        for constraint in constraints:
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

    def _locateFieldsOnNodeSetsInDofVector(self, nodeSets: list) -> dict:
        """Creates a dictionary containing the location (indices) of each entity (elements, constraints)
        within the DofVector structure.

        Returns
        -------
        dict
            A dictionary containing the location mapping.
        """

        nodeSets = nodeSets
        nodeSetFieldsInDofVector = {}

        for field in self.idcsOfNodeFieldsInDofVector:
            nodeSetFieldsInDofVector[field] = dict()
            for nSet in nodeSets:
                nodeSetFieldsInDofVector[field][nSet] = np.array(
                    [self.idcsOfFieldVariablesInDofVector[node.fields[field]] for node in nSet if field in node.fields],
                    dtype=int,
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

        sizeVIJ = self._sizeVIJ

        I = np.zeros(sizeVIJ, dtype=int)  # noqa: E741
        J = np.zeros(sizeVIJ, dtype=int)  # noqa: E741
        idxInVIJ = 0

        for (
            entity,
            entityIdcsInDofVector,
        ) in self.idcsOfHigherOrderEntitiesInDofVector.items():
            entitiesInVIJ[entity] = idxInVIJ

            nDofEntity = len(entityIdcsInDofVector)

            # looks like black magic, but it's an efficient way to generate all indices of Ke in K:
            VIJLocations = np.tile(entityIdcsInDofVector, (nDofEntity, 1))
            I[idxInVIJ : idxInVIJ + nDofEntity**2] = VIJLocations.flatten()
            J[idxInVIJ : idxInVIJ + nDofEntity**2] = VIJLocations.flatten("F")
            idxInVIJ += nDofEntity**2

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
        I = self.I  # noqa: E741
        J = self.J

        return VIJSystemMatrix(nDof, I, J, self.idcsOfHigherOrderEntitiesInVIJ)

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

        return DofVector(self.nDof, self.idcsOfHigherOrderEntitiesInDofVector)
        # return DofVector(self.nDof, self.idcsInDofVector)

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

        return self.indexToHostObjectMapping[index]

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

        if resultName not in nodeField:
            nodeField.createFieldValueEntry(resultName)

        nodeField[resultName][:] = dofVector[self.idcsOfNodeFieldsInDofVector[nodeField.name]].reshape(
            (-1, nodeField.dimension)
        )

        return nodeField

    def writeNodeFieldToDofVector(
        self, dofVector: DofVector, nodeField: NodeField, resultName: str, nodeSet: NodeSet = None
    ):
        """Write the current values of an entire NodeField to the respective locations in a given DofVector.


        Parameters
        ----------
        dofVector
            The result DofVector.
        nodeField
            The NodeField holding the values.
        resultname
            The name of the value entries held by the NodeField.
        NodeSet
            The NodeSet to consider. If None, all nodes of the NodeField are considered.

        Returns
        -------
        DofVector
            The DofVector.
        """

        if nodeSet is not None:
            if nodeSet not in self.idcsOfFieldsOnNodeSetsInDofVector[nodeField.name]:
                self.idcsOfFieldsOnNodeSetsInDofVector[nodeField.name] |= self._locateFieldsOnNodeSetsInDofVector(
                    [nodeSet]
                )
            dofVector[self.idcsOfFieldsOnNodeSetsInDofVector[nodeField.name][nodeSet]] = nodeField.subset(nodeSet)[
                resultName
            ].flatten()

        else:
            dofVector[self.idcsOfNodeFieldsInDofVector[nodeField.name]] = nodeField[resultName].flatten()

        return dofVector
