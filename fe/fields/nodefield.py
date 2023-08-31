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
from fe.config.phenomena import phenomena, getFieldSize

from fe.points.node import Node
from fe.sets.nodeset import NodeSet


class _NodeFieldValues(np.ndarray):
    """
    This class represents a result entry of NodeField,
    carrying values for each node within the field.
    It is basically a np.ndarray, but provides addional easy access by Node or NodeSets.
    The [] operator allows to access (non-contigouos read, write) at each entities location:

    Parameters
    ----------
    numberOfNodes
        The number of nodes.
    dim
        The dimension of the field.
    indicesOfNodesInArray
        The row number of each Node in the array.
    indicesOfNodeSetsInArray
        The arrays of row numbers for each NodeSet (by name) in the array.
    """

    def __new__(cls, numberOfNodes: int, dim: int, indicesOfNodesInArray: dict, indicesOfNodeSetsInArray: dict):
        obj = np.zeros((numberOfNodes, dim), dtype=float).view(cls)
        obj._indicesOfNodesInArray = indicesOfNodesInArray
        obj._indicesOfNodeSetsInArray = indicesOfNodeSetsInArray

        return obj

    def __getitem__(self, key):
        if isinstance(key, Node):
            return super().__getitem__(self._indicesOfNodesInArray[key])

        if isinstance(key, NodeSet):
            if NodeSet.label != "all":
                return super().__getitem__(self._indicesOfNodeSetsInArray[key])

        return super().__getitem__(key)

    def __setitem__(self, key, val):
        if isinstance(key, Node):
            return super().__setitem__(self._indicesOfNodesInArray[key], val)

        if isinstance(key, NodeSet):
            if NodeSet.label != "all":
                return super().__setitem__(self._indicesOfNodeSetsInArray[key], val)

        return super().__setitem__(key, val)


class NodeField:

    """
    This class represents a node field.
    A node field associates every node with multiple entries (flux, effort) of field variables.
    Furthermore, for convencience, it allows to get fast access to values for individual nodes or node sets.

    .. code-block:: console

        NodeField 'Displacement'
           values: {'U' : [[0,0,0],  # node (1)
                           [1,0,0],  # node (2)
                           [0,1,0]]} # node (3)
        (1) --------- (2)
         | *             .
         |  *+---------+   ...
         |   | [0,0,0] |
        (3)  +---------+
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
    model
        The model tree.
    """

    def __init__(self, fieldName: str, dimension: int, nodes: list, model):
        self.name = fieldName

        self.nodes = [n for n in nodes if fieldName in n.fields]

        self.dimension = dimension

        self._indicesOfNodesInArray = {n: i for i, n in enumerate(self.nodes)}
        self._indicesOfNodeSetsInArray = {
            nSet.label: np.array([self._indicesOfNodesInArray[n] for n in nSet if n in self._indicesOfNodesInArray])
            for nSet in model.nodeSets.values()
        }

        self.values = dict()

    def createFieldValueEntry(self, name: str) -> _NodeFieldValues:
        """
        Add an empty entry with given name for the field, e.g, 'U' or 'P' for flux or effort entries.

        Parameters
        ----------
        name
            The name of the entry.

        Returns
        -------
            The new entry
        """
        self.values[name] = _NodeFieldValues(
            len(self.nodes), self.dimension, self._indicesOfNodesInArray, self._indicesOfNodeSetsInArray
        )

        return self.values[name]
