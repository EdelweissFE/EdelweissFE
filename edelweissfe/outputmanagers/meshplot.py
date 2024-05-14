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
#  Magdalena Schreter magdalena.schreter@uibk.ac.at
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
@author: Magdalena

Module meshplot divided into classes:
    * Triangulation:
        creates triangles out of rectangles in order to use matplotlib plotting
        options
    * Plotter:
        creates figures, axes, grid, labels
    * Outputmanager:
        creats the plotting specific for the defined keyword lines
"""
from distutils.util import strtobool

import matplotlib.tri as mtri
import numpy as np
from matplotlib import colors

from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from edelweissfe.sets.elementset import ElementSet
from edelweissfe.utils.math import createMathExpression
from edelweissfe.utils.meshtools import (
    extractNodeCoordinatesFromElset,
    transferElsetResultsToElset,
)

documentation = {
    "figure": "figure number, (default=1)",
    "axSpec": "axis specification according to matplotlib syntax, (default=111)",
    "create=perNode": "result per node is plotted in a meshplot",
    "create=perElement": "result per element is plotted in a meshplot",
    "create=meshOnly": "plot the mesh only",
    "create=xyData": "2D Plot of results",
}


class Triangulation:
    """class that provides the division of quadrilateral elements into triangles"""

    def __init__(self, xCoord, yCoord, elNodesIdxList):
        self.triangleIdx, self.triang = self.quadIdxToTriIdx(xCoord, yCoord, elNodesIdxList)

    def quadIdxToTriIdx(self, xCoord, yCoord, elementIdxMatrix):
        triangleIdx = np.asarray(
            [[list(elIdx[:3]), [elIdx[2], elIdx[3], elIdx[0]]] for elIdx in elementIdxMatrix]
        ).reshape(-1, 3)
        triang = mtri.Triangulation(xCoord, yCoord, triangleIdx)
        return triangleIdx, triang

    def quadFieldToTriField(self, fieldValues):
        return np.asarray([2 * [fieldValues]]).reshape(-1, 1, order="f").flatten()


class MeshPlot:
    def __init__(self, coordinates, elNodesIdxList, elCoordinatesList):
        self.coordinates = coordinates
        self.elCoordinatesList = elCoordinatesList
        self.xCoord = coordinates[:, 0]
        self.yCoord = coordinates[:, 1]
        self.xLimits = [-self.xCoord.max() * 0.1, self.xCoord.max() * 1.1]
        self.yLimits = [-self.yCoord.max() * 0.1, self.yCoord.max() * 1.1]
        self.TriangObj = Triangulation(self.xCoord, self.yCoord, elNodesIdxList)
        self.contourPlotScaling = 50
        self.userColorMap = "coolwarm"

    def contourPlotFieldVariable(self, fieldValues, fig, ax, label):
        """divide quad elements into two triangles and apply a constant field
        value for both triangles"""
        resultPerTriElement = self.TriangObj.quadFieldToTriField(fieldValues)
        mapping = ax.tripcolor(
            self.xCoord,
            self.yCoord,
            self.TriangObj.triangleIdx,
            facecolors=resultPerTriElement,
            cmap=self.userColorMap,
            norm=colors.Normalize(vmax=np.nanmax(resultPerTriElement), vmin=np.nanmin(resultPerTriElement)),
        )
        cbar = fig.colorbar(mapping, fraction=0.046, pad=0.04)
        cbar.set_label(label)
        ax.set_xlim(self.xLimits)
        ax.set_ylim(self.yLimits)

    def contourPlotNodalValues(self, z, fig, ax, label, elements, nSet):
        """divide quads into two triangles and apply a nodal value to the corner nodes"""

        # PERFORMANCE IMPROVEMENT PENDING

        nSetList = [node.label for node in nSet]
        coordinates = np.asarray([node.coordinates for node in nSet])

        triangleNodes = []
        counter = 0
        for element in elements.values():
            nodeList = []
            for node in element.nodes:
                nodeList.append(node)

            if set([node.label for node in nodeList]) < set(nSetList):
                counter += 1
                triangleNodes.append([nSetList.index(node.label) for node in nodeList])

        triangObjTemp = Triangulation(coordinates[:, 0], coordinates[:, 1], triangleNodes)
        mapping = ax.tricontourf(
            triangObjTemp.triang,
            z,
            self.contourPlotScaling,
            cmap=self.userColorMap,
            norm=colors.Normalize(vmax=np.nanmax(z), vmin=np.nanmin(z)),
        )
        cbar = fig.colorbar(mapping, fraction=0.046, pad=0.04)
        cbar.set_label(label)
        ax.set_xlim(self.xLimits)
        ax.set_ylim(self.yLimits)

    def plotNodeLabels(self, labels, ax):
        """label nodes of elements"""
        for label in labels:
            ax.annotate(
                "%i" % label,
                xy=self.coordinates[label - 1, :],
                fontsize=6,
                textcoords="data",
            )

    def plotElementLabels(self, ax, elementList):
        """label nodes of elements"""
        for element in elementList:
            xCenter = 0
            yCenter = 0
            for node in element.nodes:
                xCenter += node.coordinates[0]
                yCenter += node.coordinates[1]
            ax.annotate(
                "%i" % element.elNumber,
                xy=[xCenter / 4, yCenter / 4],
                fontsize=6,
                textcoords="data",
            )

    def plotMeshGrid(self, ax, coordinateList):
        """plot grid of elements; so far only implemented for quads"""
        for element in coordinateList:
            ax.plot(
                np.append(element[:, 0], element[0, 0]),
                np.append(element[:, 1], element[0, 1]),
                "k",
                linewidth=0.3,
            )
        ax.set_xlim(self.xLimits)
        ax.set_ylim(self.yLimits)
        ax.grid(False)


class OutputManager(OutputManagerBase):
    identification = "meshPlot"

    def __init__(self, name, model, fieldOutputController, journal, plotter):
        self.domainSize = model.domainSize
        self.plotter = plotter
        self.journal = journal

        self.nodes = model.nodes
        self.elements = model.elements
        self.elSets = model.elementSets
        self.nSets = model.nodeSets

        # write List of nodeLabels
        self.labelList = np.asarray([nodeNumber for nodeNumber in self.nodes.keys()])
        # write List of node coordiantes
        self.coordinateList = np.asarray([node.coordinates for node in self.nodes.values()])
        # write list of element coordinates with 4x2 arrays (xCol, yCol)
        self.elCoordinatesList = []
        # write list of node indices for each element relevant for the meshplot output
        # in case of an 8-node element only the 4 first nodes are relevant
        self.elNodesIdxList = []

        self.perNodeJobs = []
        self.perElementJobs = []
        #        self.configJobs = []
        self.xyJobs = []
        self.saveJobs = []
        self.meshOnlyJobs = []

        self.fieldOutputController = fieldOutputController

        # # needed for meshOnly plot

    def updateDefinition(self, **kwargs: dict):
        fieldOutputController = self.fieldOutputController
        if "saveFigure" in kwargs:
            saveJob = {}
            saveJob["figure"] = kwargs.get("figure", "1")
            saveJob["fileName"] = kwargs.get("name", "exportFigure")
            saveJob["width"] = kwargs.get("width", 469.47)
            saveJob["scale"] = kwargs.get("scale", 1.0)
            saveJob["heightRatio"] = kwargs.get("heightRatio", False)
            saveJob["png"] = kwargs.get("png", True)
            self.saveJobs.append(saveJob)

        if "create" in kwargs:
            varType = kwargs["create"]

            if varType == "perNode":
                perNodeJob = {}
                perNodeJob["fieldOutput"] = kwargs["fieldOutput"]

                if type(perNodeJob["fieldOutput"].associatedSet) is ElementSet:
                    nSet = perNodeJob["fieldOutput"].associatedSet.extractNodeSet()
                else:
                    raise Exception("perNode job must be defined on a perElement fieldOutput")

                perNodeJob["nSet"] = nSet
                perNodeJob["label"] = kwargs.get("label", kwargs["fieldOutput"])
                perNodeJob["axSpec"] = kwargs.get("axSpec", "111")
                perNodeJob["figure"] = kwargs.get("figure", "1")
                perNodeJob["plotMeshGrid"] = kwargs.get("plotMeshGrid", "undeformed")
                if "f(x)" in kwargs:
                    perNodeJob["f(x)"] = createMathExpression(kwargs["f(x)"])

                # perNodeJob["dimensions"] = edelweissfe.config.phenomena.getFieldSize(
                #     perNodeJob["fieldOutput"].field, self.domainSize
                # # )
                self.perNodeJobs.append(perNodeJob)

            elif varType == "perElement":
                perElementJob = {}
                perElementJob["label"] = kwargs.get("label", kwargs["fieldOutput"])
                perElementJob["axSpec"] = kwargs.get("axSpec", "111")
                perElementJob["figure"] = kwargs.get("figure", "1")
                perElementJob["fieldOutput"] = kwargs["fieldOutput"]
                if "f(x)" in kwargs:
                    perElementJob["f(x)"] = createMathExpression(kwargs["f(x)"])

                perElementJob["plotMeshGrid"] = kwargs.get("plotMeshGrid", "unDeformed")
                self.perElementJobs.append(perElementJob)

            elif varType == "xyData":
                xyJob = {}

                if kwargs["x"] != "time":
                    xyJob["x"] = fieldOutputController.fieldOutputs[kwargs["x"]]
                else:
                    xyJob["x"] = "time"

                xyJob["y"] = fieldOutputController.fieldOutputs[kwargs["y"]]

                if "f(x)" in kwargs:
                    xyJob["f(x)"] = createMathExpression(kwargs["f(x)"])
                if "f(y)" in kwargs:
                    xyJob["f(y)"] = createMathExpression(kwargs["f(y)"], symbol="y")

                xyJob["figure"] = kwargs.get("figure", "1")
                xyJob["label"] = kwargs.get("label", xyJob["y"].name)
                xyJob["axSpec"] = kwargs.get("axSpec", "111")
                xyJob["integral"] = strtobool(kwargs.get("integral", "False"))
                self.xyJobs.append(xyJob)

            elif varType == "meshOnly":
                meshOnlyJob = {}

                meshOnlyJob["configuration"] = kwargs.get("configuration", "undeformed")

                if meshOnlyJob["configuration"] == "deformed":
                    meshOnlyJob["warpBy"] = fieldOutputController.fieldOutputs[kwargs["warpBy"]]
                meshOnlyJob["scaleFactor"] = float(kwargs.get("scaleFactor", 1.0))
                meshOnlyJob["axSpec"] = kwargs.get("axSpec", "111")
                meshOnlyJob["figure"] = kwargs.get("figure", "1")
                meshOnlyJob["plotNodeLabels"] = kwargs.get("plotNodeLabels", False)
                meshOnlyJob["plotElementLabels"] = kwargs.get("plotElementLabels", False)
                self.meshOnlyJobs.append(meshOnlyJob)

        # Initialize instance of plotterclass

    # def initializeSimulation(self, model):
    #     pass

    def initializeJob(self):
        pass

    def initializeStep(self, step):
        if self.perElementJobs or self.perNodeJobs or self.meshOnlyJobs:
            self.elCoordinatesList = extractNodeCoordinatesFromElset(self.elements.values())
            for element in self.elements.values():
                nodeIdxArray = [nodeNumber.label - 1 for nodeNumber in element.nodes[:]][:4]
                self.elNodesIdxList.append(nodeIdxArray)

            self.meshPlot = MeshPlot(self.coordinateList, self.elNodesIdxList, self.elCoordinatesList)

    def finalizeIncrement(self, **kwargs):
        pass

    def finalizeFailedIncrement(self, **kwargs):
        pass

    def finalizeStep(
        self,
    ):
        pass

    def finalizeJob(
        self,
    ):
        for xyJob in self.xyJobs:
            y = xyJob["y"].getResultHistory()

            if xyJob["x"] == "time":
                x = xyJob["y"].getTimeHistory()
            else:
                x = xyJob["x"].getResultHistory()

            if "f(x)" in xyJob:
                x = xyJob["f(x)"](x)
            if "f(y)" in xyJob:
                y = xyJob["f(y)"](y)

            self.plotter.plotXYData(x, y, xyJob["figure"], xyJob["axSpec"], xyJob)
            ax = self.plotter.getAx(xyJob["figure"], xyJob["axSpec"])
            if xyJob["integral"]:
                integral = np.trapz(y.flatten(), x=x.flatten())
                ax.fill_between(x.flatten(), 0, y.flatten(), color="gray", label=str(integral))

        for perNodeJob in self.perNodeJobs:
            result = perNodeJob["fieldOutput"].getLastResult()
            fig = self.plotter.getFig(perNodeJob["figure"])
            ax = self.plotter.getAx(perNodeJob["figure"], perNodeJob["axSpec"])
            ax.set_axis_off()
            ax.set_aspect("equal")

            if "f(x)" in perNodeJob:
                result = perNodeJob["f(x)"](result)

            if perNodeJob["plotMeshGrid"] == "undeformed":
                self.meshPlot.plotMeshGrid(ax, self.elCoordinatesList)

            result = np.squeeze(result)

            self.meshPlot.contourPlotNodalValues(
                result, fig, ax, perNodeJob["label"], self.elements, perNodeJob["nSet"]
            )

        for perElementJob in self.perElementJobs:
            fig = self.plotter.getFig(perElementJob["figure"])
            ax = self.plotter.getAx(perElementJob["figure"], perElementJob["axSpec"])
            ax.set_axis_off()
            ax.set_aspect("equal")

            #            print(self.elCoordinatesList)
            #            if perElementJob['configuration'] == 'deformed':
            #                elCoordinatesListDeformed = extractNodeCoordinatesFromElset(self.elSets['all'], perElementJob['fieldOutput'].getLastResult(), perElementJob['scaleFactor'])
            #                self.meshPlot.plotMeshGrid( ax, elCoordinatesListDeformed)

            #            self.meshPlot.plotMeshGrid(ax,  self.elCoordinatesList)

            resultArray = perElementJob["fieldOutput"].getLastResult()

            if "f(x)" in perElementJob:
                resultArray = perElementJob["f(x)"](resultArray)

            if perElementJob["fieldOutput"].associatedSet.name != "all":
                shape = (
                    (len(self.elSets["all"]), resultArray.shape[-1])
                    if resultArray.ndim >= 2
                    else len(self.elSets["all"])
                )
                resultsTarget = np.empty(shape)
                resultsTarget[:] = np.nan
                transferElsetResultsToElset(
                    self.elSets["all"],
                    perElementJob["fieldOutput"].elSet,
                    resultsTarget,
                    resultArray,
                )
                resultArray = resultsTarget

            self.meshPlot.contourPlotFieldVariable(resultArray, fig, ax, perElementJob["label"])
        for meshOnlyJob in self.meshOnlyJobs:
            fig = self.plotter.getFig(meshOnlyJob["figure"])
            ax = self.plotter.getAx(meshOnlyJob["figure"], meshOnlyJob["axSpec"])
            ax.set_axis_off()
            ax.set_aspect("equal")

            if meshOnlyJob["plotNodeLabels"]:
                self.meshPlot.plotNodeLabels(self.nodes.keys(), ax)

            if meshOnlyJob["plotElementLabels"]:
                self.meshPlot.plotElementLabels(ax, self.elSets["all"])

            if meshOnlyJob["configuration"] == "deformed":
                elCoordinatesListDeformed = extractNodeCoordinatesFromElset(
                    self.elSets["all"],
                    meshOnlyJob["warpBy"].getLastResult(),
                    meshOnlyJob["scaleFactor"],
                )
                self.meshPlot.plotMeshGrid(ax, elCoordinatesListDeformed)

            else:
                elCoordinatesListUnDeformed = extractNodeCoordinatesFromElset(self.elSets["all"])
                self.meshPlot.plotMeshGrid(ax, elCoordinatesListUnDeformed)

        #        for configJob in self.configJobs:
        #            self.plotter.configAxes(**configJob)

        for saveJob in self.saveJobs:
            self.plotter.exportFigure(
                saveJob["fileName"],
                saveJob["figure"],
                saveJob["width"],
                saveJob["scale"],
                saveJob["heightRatio"],
                saveJob["png"],
            )


##
