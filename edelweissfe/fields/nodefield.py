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

import numpy as np

from edelweissfe.points.node import Node
from edelweissfe.sets.elementset import ElementSet
from edelweissfe.sets.nodeset import NodeSet


class NodeFieldSubset:
    pass


class NodeField:
    """
    This class represents a node field.
    A node field associates every node with multiple entries (e.g., flux and effort) of a field variable.
    Furthermore, for convencience, it allows to get fast access to values for individual nodes or node sets.

    The purpose is to store field data in an efficient, contiguos manner rather than distributing it across
    all individual nodes.

    .. code-block:: console

        Example:

        NodeField 'Displacement'
           values: {'U' : [[0,0,0],  # Node (1)
                           [1,0,0],  # Node (2)
                           [0,1,0]]} # Node (3)

        Spatial representation:

        (1) --------- (2)
         | *             *
         |  *+---------+   *+---------+
         |   | [0,0,0] |    | [1,0,0] |
        (3)  +---------+    +---------+
          *
           * +---------+
             | [0,1,0] |
             +---------+

    Parameters
    ----------
    fieldName
        The name of the field.
    dimension
        The dimension of the field.
    nodes
        The associated nodes. Only nodes with active fields are considered.
    """

    def __init__(self, fieldName: str, dimension: int, nodeSet: NodeSet):
        self.name = fieldName
        self.associatedSet = nodeSet
        self.nodes = [n for n in nodeSet if fieldName in n.fields]
        self.dimension = dimension
        self._indicesOfNodesInArray = {n: i for i, n in enumerate(self.nodes)}
        self._subsetCache = dict()
        self._values = dict()

    def _getNodeFieldSubsetClass(
        self,
    ):
        return NodeFieldSubset

    def __getitem__(self, key):
        return self._values[key]

    def __contains__(self, key):
        return key in self._values

    def createFieldValueEntry(self, name: str) -> np.ndarray:
        """
        Add an empty entry with given name for the field, e.g, 'U' or 'P' for flux or effort entries.

        Parameters
        ----------
        name
            The name of the entry.

        Returns
        -------
        np.ndarray
            The new entry
        """
        self._values[name] = np.zeros((len(self.nodes), self.dimension), dtype=float)

        return self[name]

    def subset(self, subset) -> NodeFieldSubset:
        """
        Get a view on a subset of the field.

        Parameters
        ----------
        subset
            The subset, e.g., a single :class:`Node, or a :class:`NodeSet or :class:`ElementSet.

        Returns
        -------
        NodeFieldSubset
            The subset of the present NodeField.
        """
        return self._getSubsetFromCache(subset)

    def _getSubsetFromCache(self, subset) -> NodeFieldSubset:
        """
        Exploit a cache to reuse already constructed NodeFieldSubsets.
        If the subset does not exist, it will be created here.

        Parameters
        ----------
        subset
            The subset, e.g., a single Node, or a NodeSet or ElementSet.

        Returns
        -------
        NodeFieldSubset
            The subset of the present NodeField.
        """

        if subset in self._subsetCache:
            return self._subsetCache[subset]
        else:
            self._subsetCache[subset] = self._getNodeFieldSubsetClass()(self, subset)
            return self._subsetCache[subset]

    def copyEntriesFromOther(self, other, fieldValueEntries: list[str] = None):
        """
        Copy values from another NodeField.
        If the fields differ, the intersection is considered.

        Parameters
        ----------
        subset
            The sub NodeField.
        fieldValueEntries
            The list of entries which should be copied. Default: all entries are copied.
        """

        if not fieldValueEntries:
            fieldValueEntries = self._values.keys() & other._values.keys()

        commonNodes = self._indicesOfNodesInArray.keys() & other._indicesOfNodesInArray.keys()

        for fieldValueEntry in fieldValueEntries:
            self[fieldValueEntry][:] = 0.0
            idcsHere = [self._indicesOfNodesInArray[n] for n in commonNodes]
            idcsOther = [other._indicesOfNodesInArray[n] for n in commonNodes]
            self[fieldValueEntry][idcsHere] = other[fieldValueEntry][idcsOther]

    def addEntriesFromOther(self, other, fieldValueEntries: list[str] = None):
        """
        Add values from another NodeField.
        If the fields differ, the intersection is considered.

        Parameters
        ----------
        subset
            The sub NodeField.
        fieldValueEntries
            The list of entries which should be copied. Default: all entries are copied.
        """

        if not fieldValueEntries:
            fieldValueEntries = self._values.keys() & other._values.keys()

        commonNodes = self._indicesOfNodesInArray.keys() & other._indicesOfNodesInArray.keys()

        for fieldValueEntry in fieldValueEntries:
            idcsHere = [self._indicesOfNodesInArray[n] for n in commonNodes]
            idcsOther = [other._indicesOfNodesInArray[n] for n in commonNodes]
            self[fieldValueEntry][idcsHere] += other[fieldValueEntry][idcsOther]


class NodeFieldSubset(NodeField):
    def __init__(self, parentNodeField, subset):
        self.parentNodeField = parentNodeField
        self.associatedSet = subset
        self.nodes = self._getSubsetNodes(subset)
        self._indicesOfNodesInParentArray = np.array([parentNodeField._indicesOfNodesInArray[n] for n in self.nodes])

    def __getitem__(self, key):
        return self.parentNodeField[key][self._indicesOfNodesInParentArray]

    def __contains__(self, key):
        return key in self.subsetNodes

    def createFieldValueEntry(self, name):
        raise Exception("Invalid operation on subset of a NodeField")

    def subset(self, subset):
        raise Exception("Subsets of subsets are not yet implemented!")

    def _getSubsetNodes(self, subset) -> list[Node]:
        """
        Get the nodes associated with a subset.
        Only nodes with the active field are considered.

        Parameters
        ----------
        subset
            The subset, e.g., a single Node, a NodeSet or ElementSet.

        Returns
        -------
        list[Node]
            The list of subset nodes.
        """
        if type(subset) is Node:
            nodeCandidates = [
                subset,
            ]
        elif type(subset) is ElementSet:
            nodeCandidates = subset.extractNodeSet()
        elif type(subset) is NodeSet:
            nodeCandidates = subset
        else:
            raise Exception("Invalid subset")

        return [n for n in nodeCandidates if n in self.parentNodeField._indicesOfNodesInArray]
