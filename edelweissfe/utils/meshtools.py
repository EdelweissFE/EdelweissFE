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
"""
Created on Sun Jul 23 21:03:23 2017

@author: Matthias Neuner
"""
from collections import defaultdict

import numpy as np

from edelweissfe.sets.nodeset import NodeSet


def extractNodesFromElementSet(elementSet):
    """
    extract all nodes (without duplicates) from an elSet
    """
    nodeCounter = 0
    partNodes = dict()  # node -> index in nodelist
    for element in elementSet:
        for node in element.nodes:
            # if the node is already in the dict, get its index,
            # else insert it, and get the current idx = counter. increase the counter
            idx = partNodes.setdefault(node, nodeCounter)
            if idx == nodeCounter:
                # the node was just inserted, so increase the counter of inserted nodes
                nodeCounter += 1
    return NodeSet(elementSet.name, partNodes.keys())


def disassembleElsetToEnsightShapes(elementSet):
    """
    elset -> {shape : [element-index in elset, ... ], }
    """
    elements = defaultdict(list)
    for i, el in enumerate(elementSet):
        elements[el.ensightType].append(i)
    return elements


def transferElsetResultsToElset(elsetTarget, elsetOrigin, resultsTarget, resultsOrigin):
    """
    Copy results from a (sub) elSet to another elSet.
    ATTENTION: All elements in the origin set must be present in the target set. ( -> can be improved in future)
    """

    for i, el in enumerate(elsetTarget):
        el.__index__in__elsetTarget = i
    indices = [el.__index__in__elsetTarget for el in elsetOrigin]

    if resultsOrigin.ndim == 1:
        resultsTarget[indices] = resultsOrigin
    elif resultsOrigin.ndim == 2:
        resultsTarget[indices, :] = resultsOrigin

    for el in elsetTarget:
        del el.__index__in__elsetTarget


def extractNodeCoordinatesFromElset(elementSet, displacementResult=False, displacementScaleFactor=1.0, numberOfNodes=4):
    """write (deformed or undeformed) coordinates of elementSet in list format:
    [ [x1 y1 x2 y2 ... xNumberOfNodes, yNumberOfNodes], # element 1 in elementSet
      [x1 y1 x2 y2 ... xNumberOfNodes, yNumberOfNodes], # element 2 in elementSet
      ....
    ]
    """
    elCoordinatesList = []

    for element in elementSet:

        # TODO: get rid of try-except block
        try:
            nodeArray = [
                node.coordinates + displacementResult[node.label - 1, :] * displacementScaleFactor
                for node in element.nodes
            ][:numberOfNodes]
        except Exception:
            nodeArray = [node.coordinates for node in element.nodes][:numberOfNodes]
        elCoordinatesList.append(np.asarray(nodeArray))

    return elCoordinatesList
