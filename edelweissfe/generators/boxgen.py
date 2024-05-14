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
#  Paul Hofer paul.hofer@uibk.ac.at
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
A mesh generator for cuboid geometries and structured hex meshes:

.. code-block:: console

                          __ __ __ __
                        /__/__/__/__/|    A                               back
                       /__/__/__/__/ |    |                           top /
                      /__/__/__/__/ /|    | lY                         | /
                      |__|__|__|__|//|    | nY elements                |/
           y          |__|__|__|__|//|    |                    left----/----right
           |          |__|__|__|__|//|    V                           /|
           |___x      |__|__|__|__|//|   /                           / |
          /           |__|__|__|__|//   / lZ                        / bottom
         z            |__|__|__|__|/   /  nZ elements           front

                     <----lX----->
                      nX elements

nSets, elSets, surface : 'name'_left, _right, _bottom, _top, _front, _back, _all,
elSets : 'name'_centralFrontToBack, _shearBandFrontToBack, _shearBandCenterFrontToBack
are automatically generated

.. code-block:: edelweiss
    :caption: Generate meshes on the fly. Example:

    *job, name=job, domain=3d, solver=NIST

    *modelGenerator, generator=boxGen, name=gen
        nX      =4
        nY      =8
        nZ      =2
        lX      =20
        lY      =40
        lZ      =1
        elType  =C3D20R
"""

from operator import attrgetter

import numpy as np

from edelweissfe.config.elementlibrary import getElementClass
from edelweissfe.models.femodel import FEModel
from edelweissfe.points.node import Node
from edelweissfe.sets.elementset import ElementSet
from edelweissfe.sets.nodeset import NodeSet
from edelweissfe.utils.misc import convertLinesToStringDictionary

documentation = {
    "x0": "(optional) origin at x axis",
    "y0": "(optional) origin at y axis",
    "z0": "(optional) origin at z axis",
    "lX": "(optional) length of the body along x axis",
    "lY": "(optional) length of the body along y axis",
    "lZ": "(optional) length of the body along z axis",
    "nX": "(optional) number of elements along x",
    "nY": "(optional) number of elements along y",
    "nZ": "(optional) number of elements along z",
    "elType": "type of element",
}


def generateModelData(generatorDefinition: dict, model: FEModel, journal) -> dict:
    options = generatorDefinition["data"]
    options = convertLinesToStringDictionary(options)

    name = generatorDefinition.get("name", "boxGen")

    x0 = float(options.get("x0", 0.0))
    y0 = float(options.get("y0", 0.0))
    z0 = float(options.get("z0", 0.0))
    lX = float(options.get("lX", 1.0))
    lY = float(options.get("lY", 1.0))
    lZ = float(options.get("lZ", 1.0))
    nX = int(options.get("nX", 1))
    nY = int(options.get("nY", 1))
    nZ = int(options.get("nZ", 1))
    elType = getElementClass(options.get("elProvider", None))

    testEl = elType(
        options["elType"],
        0,
    )

    if testEl.nNodes == 8:
        nNodesX = nX + 1
        nNodesY = nY + 1
        nNodesZ = nZ + 1
    elif testEl.nNodes == 20:
        nNodesX = 2 * nX + 1
        nNodesY = 2 * nY + 1
        nNodesZ = 2 * nZ + 1
    else:
        return

    # coordinates of layers
    xLayers = np.linspace(x0, x0 + lX, nNodesX)
    yLayers = np.linspace(y0, y0 + lY, nNodesY)
    zLayers = np.linspace(z0, z0 + lZ, nNodesZ)

    nodes = []
    currentNodeLabel = 1
    if model.nodes:
        currentNodeLabel += max(model.nodes.keys())
    for ix in range(nNodesX):
        for iy in range(nNodesY):
            for iz in range(nNodesZ):
                node = Node(currentNodeLabel, np.array([xLayers[ix], yLayers[iy], zLayers[iz]]))
                nodes.append(node)
                # only add node to model if it will be part of an element
                if testEl.nNodes == 8 or testEl.nNodes == 20 and sum(np.mod([ix, iy, iz], 2)) < 2:
                    model.nodes[currentNodeLabel] = node
                    currentNodeLabel += 1

    # # 3d plot of nodes; for debugging
    # def plotNodeList( nodeList ):
    #     nodeListFile = "nodes.dat"
    #     with open( nodeListFile, "w+" ) as f:
    #         for node in nodeList:
    #             coords = node.coordinates
    #             line = "{:5}, {:12}, {:12}, {:12}\n".format( node.label, coords[0], coords[1], coords[2] )
    #             f.write( line )

    #     cmd = [ "gnuplot",
    #             "plotConfig",
    #             "-p",
    #             "-e",
    #             "\' filename=\"{}\"; splot filename using 4:2:3:(sprintf(\"(%i)\", $1)) with labels \'".format( nodeListFile ) ]
    #     os.system( " ".join( cmd ) )

    # plotNodeList( nodes )
    # plotNodeList( [model.nodes[n] for n in model.nodes] )

    # fmt: off

    elements = []
    currentElementLabel = 1
    if model.elements:
        currentElementLabel += max(model.elements.keys())
    for ix in range(nX):
        for iy in range(nY):
            for iz in range(nZ):
                if testEl.nNodes == 8:
                    nodeList = [
                        nodes[0 + ix * (nNodesY * nNodesZ) + iy * nNodesZ + iz],
                        nodes[1 + ix * (nNodesY * nNodesZ) + iy * nNodesZ + iz],
                        nodes[1 + (1 + ix) * (nNodesY * nNodesZ) + iy * nNodesZ + iz],
                        nodes[0 + (1 + ix) * (nNodesY * nNodesZ) + iy * nNodesZ + iz],
                        nodes[0 + ix * (nNodesY * nNodesZ) + (1 + iy) * nNodesZ + iz],
                        nodes[1 + ix * (nNodesY * nNodesZ) + (1 + iy) * nNodesZ + iz],
                        nodes[
                            1 + (1 + ix) * (nNodesY * nNodesZ) + (1 + iy) * nNodesZ + iz
                        ],
                        nodes[
                            0 + (1 + ix) * (nNodesY * nNodesZ) + (1 + iy) * nNodesZ + iz
                        ],
                    ]
                elif testEl.nNodes == 20:
                    nodeList = [
                        nodes[
                            0 + 2 * ix * (nNodesY * nNodesZ) + 2 * iy * nNodesZ + 2 * iz
                        ],
                        nodes[
                            2 + 2 * ix * (nNodesY * nNodesZ) + 2 * iy * nNodesZ + 2 * iz
                        ],
                        nodes[
                            2
                            + 2 * (1 + ix) * (nNodesY * nNodesZ)
                            + 2 * iy * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            0
                            + 2 * (1 + ix) * (nNodesY * nNodesZ)
                            + 2 * iy * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            0
                            + 2 * ix * (nNodesY * nNodesZ)
                            + 2 * (1 + iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            2
                            + 2 * ix * (nNodesY * nNodesZ)
                            + 2 * (1 + iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            2
                            + 2 * (1 + ix) * (nNodesY * nNodesZ)
                            + 2 * (1 + iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            0
                            + 2 * (1 + ix) * (nNodesY * nNodesZ)
                            + 2 * (1 + iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            1 + 2 * ix * (nNodesY * nNodesZ) + 2 * iy * nNodesZ + 2 * iz
                        ],
                        nodes[
                            2
                            + (1 + 2 * ix) * (nNodesY * nNodesZ)
                            + 2 * iy * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            1
                            + 2 * (1 + ix) * (nNodesY * nNodesZ)
                            + 2 * iy * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            0
                            + (1 + 2 * ix) * (nNodesY * nNodesZ)
                            + 2 * iy * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            1
                            + 2 * ix * (nNodesY * nNodesZ)
                            + 2 * (1 + iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            2
                            + (1 + 2 * ix) * (nNodesY * nNodesZ)
                            + 2 * (1 + iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            1
                            + 2 * (1 + ix) * (nNodesY * nNodesZ)
                            + 2 * (1 + iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            0
                            + (1 + 2 * ix) * (nNodesY * nNodesZ)
                            + 2 * (1 + iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            0
                            + 2 * ix * (nNodesY * nNodesZ)
                            + (1 + 2 * iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            2
                            + 2 * ix * (nNodesY * nNodesZ)
                            + (1 + 2 * iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            2
                            + 2 * (1 + ix) * (nNodesY * nNodesZ)
                            + (1 + 2 * iy) * nNodesZ
                            + 2 * iz
                        ],
                        nodes[
                            0
                            + 2 * (1 + ix) * (nNodesY * nNodesZ)
                            + (1 + 2 * iy) * nNodesZ
                            + 2 * iz
                        ],
                    ]
                else:
                    return

                # plotNodeList( nodeList )

                # newEl = elType(options["elType"], nodeList, currentElementLabel)
                newEl = elType(options["elType"], currentElementLabel)
                newEl.setNodes(nodeList)

                elements.append(newEl)
                model.elements[currentElementLabel] = newEl

                # for i, node in enumerate(newEl.nodes):
                #     node.fields.update([(f, True) for f in newEl.fields[i]])

                currentElementLabel += 1

    # fmt: on
    # model.initializeNodeFields()
    model._populateNodeFieldVariablesFromElements()

    nG = np.asarray(nodes).reshape(nNodesX, nNodesY, nNodesZ)

    # nodesets:
    nodeSets = []

    # 6 faces
    getFields = np.vectorize(attrgetter("fields"))
    getLength = np.vectorize(len)
    filterGrid = getLength(getFields(nG)) > 0

    def getFilteredNodes(s):
        return nG[s][filterGrid[s]]

    nodeSets.append(NodeSet("{:}_top".format(name), getFilteredNodes(np.s_[:, -1, :])))
    nodeSets.append(NodeSet("{:}_bottom".format(name), getFilteredNodes(np.s_[:, 0, :])))
    nodeSets.append(NodeSet("{:}_right".format(name), getFilteredNodes(np.s_[-1, :, :])))
    nodeSets.append(NodeSet("{:}_left".format(name), getFilteredNodes(np.s_[0, :, :])))
    nodeSets.append(NodeSet("{:}_front".format(name), getFilteredNodes(np.s_[:, :, -1])))
    nodeSets.append(NodeSet("{:}_back".format(name), getFilteredNodes(np.s_[:, :, 0])))

    # 12 edges
    nodeSets.append(NodeSet("{:}_bottomRight".format(name), getFilteredNodes(np.s_[-1, 0, :])))
    nodeSets.append(NodeSet("{:}_bottomLeft".format(name), getFilteredNodes(np.s_[0, 0, :])))
    nodeSets.append(NodeSet("{:}_bottomFront".format(name), getFilteredNodes(np.s_[:, 0, -1])))
    nodeSets.append(NodeSet("{:}_bottomBack".format(name), getFilteredNodes(np.s_[:, 0, 0])))

    nodeSets.append(NodeSet("{:}_topRight".format(name), getFilteredNodes(np.s_[-1, -1, :])))
    nodeSets.append(NodeSet("{:}_topLeft".format(name), getFilteredNodes(np.s_[0, -1, :])))
    nodeSets.append(NodeSet("{:}_topFront".format(name), getFilteredNodes(np.s_[:, -1, -1])))
    nodeSets.append(NodeSet("{:}_topBack".format(name), getFilteredNodes(np.s_[:, -1, 0])))

    nodeSets.append(NodeSet("{:}_rightBack".format(name), getFilteredNodes(np.s_[-1, :, 0])))
    nodeSets.append(NodeSet("{:}_rightFront".format(name), getFilteredNodes(np.s_[-1, :, -1])))

    nodeSets.append(NodeSet("{:}_leftBack".format(name), getFilteredNodes(np.s_[0, :, 0])))
    nodeSets.append(NodeSet("{:}_leftFront".format(name), getFilteredNodes(np.s_[0, :, -1])))

    nodeSets.append(NodeSet("{:}_centerX".format(name), getFilteredNodes(np.s_[int(nNodesX / 2), :, :])))
    nodeSets.append(NodeSet("{:}_centerY".format(name), getFilteredNodes(np.s_[:, int(nNodesY / 2), :])))
    nodeSets.append(NodeSet("{:}_centerZ".format(name), getFilteredNodes(np.s_[:, :, int(nNodesZ / 2)])))

    # 8 vertices
    nodeSets.append(NodeSet("{:}_bottomRightFront".format(name), nG[-1, 0, -1]))
    nodeSets.append(NodeSet("{:}_bottomRightBack".format(name), nG[-1, 0, 0]))
    nodeSets.append(NodeSet("{:}_bottomLeftFront".format(name), nG[0, 0, -1]))
    nodeSets.append(NodeSet("{:}_bottomLeftBack".format(name), nG[0, 0, 0]))

    nodeSets.append(NodeSet("{:}_topRightFront".format(name), nG[-1, -1, -1]))
    nodeSets.append(NodeSet("{:}_topRightBack".format(name), nG[-1, -1, 0]))
    nodeSets.append(NodeSet("{:}_topLeftFront".format(name), nG[0, -1, -1]))
    nodeSets.append(NodeSet("{:}_topLeftBack".format(name), nG[0, -1, 0]))

    for nodeSet in nodeSets:
        model.nodeSets[nodeSet.name] = nodeSet

    # element sets
    elementSets = []
    elementSets.append(ElementSet("{:}_all".format(name), elements))

    elGrid = np.asarray(elements).reshape(nX, nY, nZ)
    elementSets.append(ElementSet("{:}_bottom".format(name), np.ravel(elGrid[:, 0, :])))
    elementSets.append(ElementSet("{:}_top".format(name), np.ravel(elGrid[:, -1, :])))
    elementSets.append(ElementSet("{:}_right".format(name), np.ravel(elGrid[-1, :, :])))
    elementSets.append(ElementSet("{:}_left".format(name), np.ravel(elGrid[0, :, :])))
    elementSets.append(ElementSet("{:}_front".format(name), np.ravel(elGrid[:, :, -1])))
    elementSets.append(ElementSet("{:}_back".format(name), np.ravel(elGrid[:, :, 0])))

    elementSets.append(
        ElementSet(
            "{:}_centralFrontToBack".format(name),
            np.ravel(elGrid[int(nX / 2), int(nY / 2), 0:nZ]),
        )
    )

    elementSets.append(ElementSet("{:}_centerSliceX".format(name), np.ravel(elGrid[int(nX / 2), :, :])))
    elementSets.append(ElementSet("{:}_centerSliceY".format(name), np.ravel(elGrid[:, int(nY / 2), :])))
    elementSets.append(ElementSet("{:}_centerSliceZ".format(name), np.ravel(elGrid[:, :, int(nZ / 2)])))

    nShearBand = min(nX, nY)
    if nShearBand > 3:
        shearBand = []
        for i1 in range(nShearBand):
            shearBand.extend(
                np.ravel(
                    elGrid[
                        int(nX / 2 + i1 - nShearBand / 2),
                        int(nY / 2 + i1 - nShearBand / 2),
                        0:nZ,
                    ]
                )
            )
        elementSets.append(ElementSet("{:}_shearBandFrontToBack".format(name), [e for e in shearBand]))
        elementSets.append(
            ElementSet(
                "{:}_shearBandCenterFrontToBack".format(name),
                [e for e in shearBand[(int(nShearBand / 2) - 1) * nZ : (int(nShearBand / 2) + 2) * nZ]],
            )
        )

    # model.elementSets["{:}_sandwichHorizontal".format(name)] = []
    # for elList in elGrid[1:-1, :]:
    #     for e in elList:
    #         model.elementSets["{:}_sandwichHorizontal".format(name)].append(e)

    # model.elementSets["{:}_sandwichVertical".format(name)] = []
    # for elList in elGrid[:, 1:-1]:
    #     for e in elList:
    #         model.elementSets["{:}_sandwichVertical".format(name)].append(e)

    # model.elementSets["{:}_core".format(name)] = []
    # for elList in elGrid[1:-1, 1:-1]:
    #     for e in elList:
    #         model.elementSets["{:}_core".format(name)].append(e)

    for elementSet in elementSets:
        model.elementSets[elementSet.name] = elementSet

    # surfaces
    model.surfaces["{:}_bottom".format(name)] = {1: model.elementSets["{:}_bottom".format(name)]}
    model.surfaces["{:}_top".format(name)] = {2: model.elementSets["{:}_top".format(name)]}

    model.surfaces["{:}_right".format(name)] = {5: model.elementSets["{:}_right".format(name)]}
    model.surfaces["{:}_left".format(name)] = {3: model.elementSets["{:}_left".format(name)]}

    model.surfaces["{:}_front".format(name)] = {4: model.elementSets["{:}_front".format(name)]}
    model.surfaces["{:}_back".format(name)] = {6: model.elementSets["{:}_back".format(name)]}

    return model
