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
# Created on Wed Apr 12 15:41:51 2017

# @author: Matthias Neuner
"""

A mesh generator, for rectangular geometries and structured quad meshes:


.. code-block:: console

        <-----l----->
         nX elements
         __ __ __ __
        |__|__|__|__|  A
        |__|__|__|__|  |
        |__|__|__|__|  | h
        |__|__|__|__|  | nY elements
      | |__|__|__|__|  |
      | |__|__|__|__|  V
    x0|_____
      y0

nSets, elSets, surface : 'name'_top, _bottom, _left, _right, ...
are automatically generated

Datalines:
"""

import numpy as np

from edelweissfe.config.elementlibrary import getElementClass
from edelweissfe.points.node import Node
from edelweissfe.sets.elementset import ElementSet
from edelweissfe.sets.nodeset import NodeSet
from edelweissfe.utils.misc import convertLinesToStringDictionary

documentation = {
    "x0": "(optional) origin at x axis",
    "y0": "(optional) origin at y axis",
    "h": "(optional) height of the body",
    "l": "(optional) length of the body",
    "nX": "(optional) number of elements along x",
    "nY": "(optional) number of elements along y",
    "elType": "type of element",
}


def generateModelData(generatorDefinition, model, journal):
    options = generatorDefinition["data"]
    options = convertLinesToStringDictionary(options)

    name = generatorDefinition.get("name", "planeRect")

    x0 = float(options.get("x0", 0.0))
    y0 = float(options.get("y0", 0.0))
    h = float(options.get("h", 1.0))
    l = float(options.get("l", 1.0))  # noqa: E741
    nX = int(options.get("nX", 10))
    nY = int(options.get("nY", 10))
    elType = getElementClass(options["elType"], options.get("elProvider", None))

    testEl = elType(
        options["elType"],
        0,
    )
    if testEl.nNodes == 4:
        nNodesX = nX + 1
        nNodesY = nY + 1

    if testEl.nNodes == 8:
        nNodesX = 2 * nX + 1
        nNodesY = 2 * nY + 1

    grid = np.mgrid[
        x0 : x0 + l : nNodesX * 1j,
        y0 : y0 + h : nNodesY * 1j,
    ]

    nodes = []
    currentNodeLabel = 1

    for x in range(nNodesX):
        for y in range(nNodesY):
            node = Node(currentNodeLabel, grid[:, x, y])
            model.nodes[currentNodeLabel] = node
            nodes.append(node)
            currentNodeLabel += 1

    nG = np.asarray(nodes).reshape(nNodesX, nNodesY)

    currentElementLabel = 1

    elements = []
    for x in range(nX):
        for y in range(nY):
            if testEl.nNodes == 4:
                newEl = elType(options["elType"], currentElementLabel)
                newEl.setNodes([nG[x, y], nG[x + 1, y], nG[x + 1, y + 1], nG[x, y + 1]])

            elif testEl.nNodes == 8:
                newEl = elType(
                    options["elType"],
                    currentElementLabel,
                )
                newEl.setNodes(
                    [
                        nG[2 * x, 2 * y],
                        nG[2 * x + 2, 2 * y],
                        nG[2 * x + 2, 2 * y + 2],
                        nG[2 * x, 2 * y + 2],
                        nG[2 * x + 1, 2 * y],
                        nG[2 * x + 2, 2 * y + 1],
                        nG[2 * x + 1, 2 * y + 2],
                        nG[2 * x, 2 * y + 1],
                    ]
                )
            elements.append(newEl)
            model.elements[currentElementLabel] = newEl

            currentElementLabel += 1

    model._populateNodeFieldVariablesFromElements()

    # nodesets:
    model.nodeSets["{:}_all".format(name)] = NodeSet("{:}_all".format(name), [n for n in np.ravel(nG) if len(n.fields)])

    model.nodeSets["{:}_left".format(name)] = NodeSet("{:}_left".format(name), np.ravel(nG[0, :]))
    model.nodeSets["{:}_right".format(name)] = NodeSet("{:}_right".format(name), np.ravel(nG[-1, :]))
    model.nodeSets["{:}_top".format(name)] = NodeSet("{:}_top".format(name), np.ravel(nG[:, -1]))
    model.nodeSets["{:}_bottom".format(name)] = NodeSet("{:}_bottom".format(name), np.ravel(nG[:, 0]))

    model.nodeSets["{:}_leftBottom".format(name)] = NodeSet("{:}_leftBottom".format(name), nG[0, 0])
    model.nodeSets["{:}_leftTop".format(name)] = NodeSet("{:}_leftTop".format(name), nG[0, -1])
    model.nodeSets["{:}_rightBottom".format(name)] = NodeSet("{:}_rightBottom".format(name), nG[-1, 0])
    model.nodeSets["{:}_rightTop".format(name)] = NodeSet("{:}_rightTop".format(name), nG[-1, -1])

    # element sets
    elGrid = np.asarray(elements).reshape(nX, nY)
    model.elementSets["{:}_bottom".format(name)] = ElementSet("{:}_bottom".format(name), np.ravel(elGrid[:, 0]))
    model.elementSets["{:}_top".format(name)] = ElementSet("{:}_top".format(name), np.ravel(elGrid[:, -1]))
    model.elementSets["{:}_central".format(name)] = ElementSet(
        "{:}_central".format(name), elGrid[int(nX / 2), int(nY / 2)]
    )
    model.elementSets["{:}_right".format(name)] = ElementSet("{:}_right".format(name), np.ravel(elGrid[-1, :]))
    model.elementSets["{:}_left".format(name)] = ElementSet("{:}_left".format(name), np.ravel(elGrid[0, :]))

    nShearBand = min(nX, nY)
    if nShearBand > 3:
        shearBand = [
            elGrid[int(nX / 2 + i - nShearBand / 2), int(nY / 2 + i - nShearBand / 2)] for i in range(nShearBand)
        ]
        model.elementSets["{:}_shearBand".format(name)] = ElementSet(
            "{:}_shearBand".format(name), [e for e in shearBand]
        )
        model.elementSets["{:}_shearBandCenter".format(name)] = ElementSet(
            "{:}_shearBandCenter".format(name),
            [e for e in shearBand[int(nShearBand / 2) - 1 : int(nShearBand / 2) + 2]],
        )

    model.elementSets["{:}_sandwichHorizontal".format(name)] = ElementSet(
        "{:}_sandwichHorizontal".format(name), np.ravel(elGrid[1:-1, :])
    )

    model.elementSets["{:}_sandwichVertical".format(name)] = ElementSet(
        "{:}_sandwichVertical".format(name), np.ravel(elGrid[:, 1:-1])
    )

    model.elementSets["{:}_core".format(name)] = ElementSet("{:}_core".format(name), np.ravel(elGrid[1:-1, 1:-1]))

    model.elementSets["{:}_all".format(name)] = ElementSet("{:}_all".format(name), np.ravel(elGrid))
    # surfaces
    model.surfaces["{:}_bottom".format(name)] = {1: np.ravel(elGrid[:, 0])}
    model.surfaces["{:}_top".format(name)] = {3: np.ravel(elGrid[:, -1])}
    model.surfaces["{:}_right".format(name)] = {2: np.ravel(elGrid[-1, :])}
    model.surfaces["{:}_left".format(name)] = {4: np.ravel(elGrid[0, :])}

    return model
