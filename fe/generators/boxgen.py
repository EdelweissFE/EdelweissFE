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

Datalines:
"""
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

from fe.elements.node import Node
from fe.config.elementlibrary import getElementByName
from fe.utils.misc import stringDict

import numpy as np
import os


def generateModelData(generatorDefinition, modelInfo, journal):

    options = generatorDefinition["data"]
    options = stringDict([e for l in options for e in l])

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
    elType = getElementByName(options["elType"], options.get("elProvider", None))

    testEl = elType(
        options["elType"],
        [],
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
    for ix in range(nNodesX):
        for iy in range(nNodesY):
            for iz in range(nNodesZ):
                node = Node(currentNodeLabel, [xLayers[ix], yLayers[iy], zLayers[iz]])
                nodes.append(node)
                currentNodeLabel += 1
                # only add node to modelInfo if it will be part of an element
                if (
                    testEl.nNodes == 8
                    or testEl.nNodes == 20
                    and sum(np.mod([ix, iy, iz], 2)) < 2
                ):
                    modelInfo["nodes"][currentNodeLabel] = node

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
    # plotNodeList( [modelInfo["nodes"][n] for n in modelInfo["nodes"]] )

    elements = []
    currentElementLabel = 1
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

                newEl = elType(options["elType"], nodeList, currentElementLabel)

                elements.append(newEl)
                modelInfo["elements"][currentElementLabel] = newEl

                for i, node in enumerate(newEl.nodes):
                    node.fields.update([(f, True) for f in newEl.fields[i]])

                currentElementLabel += 1

    nG = np.asarray(nodes).reshape(nNodesX, nNodesY, nNodesZ)

    # nodesets:
    modelInfo["nodeSets"]["{:}_all".format(name)] = []
    for n in np.ravel(nG):
        if len(n.fields) > 0:
            modelInfo["nodeSets"]["{:}_all".format(name)].append(n)

    modelInfo["nodeSets"]["{:}_left".format(name)] = [
        n for n in np.ravel(nG[0, :, :]) if len(n.fields) > 0
    ]
    modelInfo["nodeSets"]["{:}_right".format(name)] = [
        n for n in np.ravel(nG[-1, :, :]) if len(n.fields) > 0
    ]
    modelInfo["nodeSets"]["{:}_top".format(name)] = [
        n for n in np.ravel(nG[:, -1, :]) if len(n.fields) > 0
    ]
    modelInfo["nodeSets"]["{:}_bottom".format(name)] = [
        n for n in np.ravel(nG[:, 0, :]) if len(n.fields) > 0
    ]
    modelInfo["nodeSets"]["{:}_front".format(name)] = [
        n for n in np.ravel(nG[:, :, -1]) if len(n.fields) > 0
    ]
    modelInfo["nodeSets"]["{:}_back".format(name)] = [
        n for n in np.ravel(nG[:, :, 0]) if len(n.fields) > 0
    ]

    modelInfo["nodeSets"]["{:}_bottomLeft".format(name)] = [
            n for n in np.ravel(nG[0, 0, :]) if len(n.fields) > 0
    ]
    modelInfo["nodeSets"]["{:}_bottomRight".format(name)] = [
            n for n in np.ravel(nG[-1, 0, :]) if len(n.fields) > 0
    ]
    modelInfo["nodeSets"]["{:}_bottomFront".format(name)] = [
        n for n in np.ravel(nG[:, 0, -1]) if len(n.fields) > 0
    ]
    modelInfo["nodeSets"]["{:}_bottomBack".format(name)] = [
        n for n in np.ravel(nG[:, 0, 0]) if len(n.fields) > 0
    ]

    # modelInfo["nodeSets"]["{:}_leftBottom".format(name)]  = [nG[0, 0]]
    # modelInfo["nodeSets"]["{:}_leftTop".format(name)]     = [nG[0, -1]]
    # modelInfo["nodeSets"]["{:}_rightBottom".format(name)] = [nG[-1, 0]]
    # modelInfo["nodeSets"]["{:}_rightTop".format(name)]    = [nG[-1, -1]]

    # element sets
    elGrid = np.asarray(elements).reshape(nX, nY, nZ)
    modelInfo["elementSets"]["{:}_bottom".format(name)] = np.ravel(elGrid[:, 0, :])
    modelInfo["elementSets"]["{:}_top".format(name)] = np.ravel(elGrid[:, -1, :])
    modelInfo["elementSets"]["{:}_right".format(name)] = np.ravel(elGrid[-1, :, :])
    modelInfo["elementSets"]["{:}_left".format(name)] = np.ravel(elGrid[0, :, :])
    modelInfo["elementSets"]["{:}_front".format(name)] = np.ravel(elGrid[:, :, -1])
    modelInfo["elementSets"]["{:}_back".format(name)] = np.ravel(elGrid[:, :, 0])

    modelInfo["elementSets"]["{:}_centralFrontToBack".format(name)] = np.ravel(
        elGrid[int(nX / 2), int(nY / 2), 0:nZ]
    )

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
        modelInfo["elementSets"]["{:}_shearBandFrontToBack".format(name)] = [
            e for e in shearBand
        ]
        modelInfo["elementSets"]["{:}_shearBandCenterFrontToBack".format(name)] = [
            e
            for e in shearBand[
                (int(nShearBand / 2) - 1) * nZ : (int(nShearBand / 2) + 2) * nZ
            ]
        ]

    # modelInfo["elementSets"]["{:}_sandwichHorizontal".format(name)] = []
    # for elList in elGrid[1:-1, :]:
    #     for e in elList:
    #         modelInfo["elementSets"]["{:}_sandwichHorizontal".format(name)].append(e)

    # modelInfo["elementSets"]["{:}_sandwichVertical".format(name)] = []
    # for elList in elGrid[:, 1:-1]:
    #     for e in elList:
    #         modelInfo["elementSets"]["{:}_sandwichVertical".format(name)].append(e)

    # modelInfo["elementSets"]["{:}_core".format(name)] = []
    # for elList in elGrid[1:-1, 1:-1]:
    #     for e in elList:
    #         modelInfo["elementSets"]["{:}_core".format(name)].append(e)

    # surfaces
    modelInfo["surfaces"]["{:}_bottom".format(name)] = {
        1: [e for e in np.ravel(elGrid[:, 0, :])]
    }
    modelInfo["surfaces"]["{:}_top".format(name)] = {
        2: [e for e in np.ravel(elGrid[:, -1, :])]
    }

    modelInfo["surfaces"]["{:}_right".format(name)] = {
        5: [e for e in np.ravel(elGrid[-1, :, :])]
    }
    modelInfo["surfaces"]["{:}_left".format(name)] = {
        3: [e for e in np.ravel(elGrid[0, :, :])]
    }

    modelInfo["surfaces"]["{:}_front".format(name)] = {
        4: [e for e in np.ravel(elGrid[:, :, -1])]
    }  #
    modelInfo["surfaces"]["{:}_back".format(name)] = {
        6: [e for e in np.ravel(elGrid[:, :, 0])]
    }

    return modelInfo
