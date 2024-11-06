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

from abc import ABC, abstractmethod

import numpy as np

from edelweissfe.points.node import Node


class BaseNodeCouplingEntity(ABC):
    """
    This class represents the most basic entity which couples nodes (and their fields) in a numerical model.
    In the context of finite element analysis, this entity is is the base class for finite elements.

    The entity information is used for setting up the global system of equations.
    Each entity is defined by a list of nodes it couples, and a list of fields per node, on which it actually operates.

    The reason for a common base class is that not only finite elements couple nodes, but also constraints, cells (MPM) or particles (SPH, RKPM).
    Having a common base class allows for a unified treatment of these entities in the global system of equations.
    """

    @property
    @abstractmethod
    def nNodes(self) -> int:
        """The number of nodes this entity couples."""

    @property
    @abstractmethod
    def nodes(self) -> list[Node]:
        """The list of nodes this currently entity holds."""

    @property
    @abstractmethod
    def nDof(self) -> int:
        """The total number of degrees of freedom this entity has.
        Not that this information is redundent, as it results from the list of fields on nodes,
        but we keep this redundant information for the sake of performance."""

    @property
    @abstractmethod
    def fields(self) -> list[list[str]]:
        """The list of fields per node which this entity couples."""

    @property
    @abstractmethod
    def dofIndicesPermutation(self) -> np.ndarray | None:
        """If the provides computes residual vectors and stiffness matrices not nodewise, but e.g., fieldwise,
        this permutation pattern is used to aggregate all entries in order to resemble the defined fields nodewise.

        For performance reasons, this is computed once and stored for later use."""

    @property
    @abstractmethod
    def ensightType(self) -> str:
        """The shape of the element in Ensight Gold notation.
        Valid types are:
        point
        bar2
        bar3
        tria3
        tria6
        quad4
        quad8
        tetra4
        tetra10
        pyramid5
        pyramid13
        penta6
        penta15
        hexa8
        hexa20
        nsided
        nfaced
        """

    @property
    @abstractmethod
    def visualizationNodes(self) -> list[Node]:
        """The nodes for visualization. Commonly, these are the same as the nodes of the entity.
        However, in some cases, the visualization nodes are different from the nodes of the entity, e.g., in case of mixed formulations.
        """
