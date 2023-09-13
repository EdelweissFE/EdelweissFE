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
# Created on Mon Apr 17 11:37:26 2017

# @author: Matthias Neuner
"""

Plot result for a nodeSet or an elementSet along the true geometrical distance.
Corresponds to the plot along path functionality in Abaqus.
"""
documentation = {
    "fieldOutput": "Name of the field output to be plotted (defined on a nodeSet/elementSet)",
    "figure": "(Optional), figure of the Plotter",
    "axSpec": "(Optional), axSpec (MATLAB syntax) in the figure",
    "label": "(Optional), label, standard=fieldOutputs name",
    "f(x)": "(Optional), apply math",
    "normalize": "(Optional), normalize peak to 1.0",
}

from fe.outputmanagers.base.outputmanagerbase import OutputManagerBase

from fe.utils.misc import convertLineToStringDictionary
from fe.utils.math import createMathExpression
import numpy as np


class OutputManager(OutputManagerBase):
    identification = "PathPlotter"

    def __init__(self, name, definitionLines, model, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []
        self.plotter = plotter

        for defline in definitionLines:
            entry = {}
            defDict = convertLineToStringDictionary(defline)
            entry["fieldOutput"] = fieldOutputController.fieldOutputs[defDict["fieldOutput"]]

            # compute distance(s), entity 0 is the reference entity in the 'origin'
            entry["pathDistances"] = [0.0]
            entry["nStages"] = int(defDict.get("nStages", 1))
            entry["export"] = bool(defDict.get("export", False))

            try:  # nSet?
                nodes = entry["fieldOutput"].nSet
                # 1) distances between nodes:
                distances = [
                    np.linalg.norm(nodes[i + 1].coordinates - nodes[i].coordinates) for i in range(len(nodes) - 1)
                ]
            except AttributeError:  # no, its an elSet!
                elements = entry["fieldOutput"].elSet
                # dirty computation of centroid by taking the mean (not correct, but fast)
                elCentroids = [np.asarray(el.nodeCoordinates).reshape(el.nNodes, -1).mean(axis=0) for el in elements]
                elCentroids = np.asarray(elCentroids)
                # 1) distances between elements:
                distances = [
                    np.linalg.norm(elCentroids[i + 1, :] - elCentroids[i, :]) for i in range(len(elCentroids) - 1)
                ]

            # 2) distances with respect to entity 0
            for dist in distances:
                entry["pathDistances"].append(entry["pathDistances"][-1] + dist)

            entry["f(x)"] = createMathExpression(defDict.get("f(x)", "x"))
            entry["label"] = defDict.get("label", entry["fieldOutput"].name)

            entry["figure"] = defDict.get("figure", 1)
            entry["axSpec"] = defDict.get("axSpec", 111)

            entry["normalize"] = defDict.get("normalize", False)
            self.monitorJobs.append(entry)

            if "exportPath" in defDict:
                np.savetxt(
                    defDict["exportPath"],
                    np.asarray(entry["pathDistances"]),
                )

    # def initializeSimulation(self, model):
    #     pass

    def initializeStep(self, step):
        for nJob in self.monitorJobs:
            self.plotStages = np.linspace(0, step.length, nJob["nStages"])

    def finalizeIncrement(self, increment, **kwargs):
        totalTime = increment[3] + increment[4]
        if totalTime > self.plotStages[0]:
            for nJob in self.monitorJobs:
                nJob_ = nJob.copy()
                nJob_["label"] = None
                result = nJob["fieldOutput"].getLastResult()

                result = nJob["f(x)"](result)

                result = np.squeeze(result)

                if nJob["normalize"]:
                    result /= np.max(np.abs(result))

                if nJob["export"]:
                    exportData = np.column_stack((nJob["pathDistances"], result))
                    np.savetxt(
                        nJob["label"] + "stage_" + str(nJob["nStages"] - len(self.plotStages)) + ".csv", exportData
                    )

                self.plotter.plotXYData(nJob["pathDistances"], result, nJob["figure"], nJob["axSpec"], nJob_)

            self.plotStages = np.delete(self.plotStages, 0)

    def finalizeFailedIncrement(self, **kwargs):
        pass

    def finalizeStep(
        self,
    ):
        pass

    def finalizeJob(
        self,
    ):
        for nJob in self.monitorJobs:
            result = nJob["fieldOutput"].getLastResult()
            result = nJob["f(x)"](result)

            result = np.squeeze(result)

            if nJob["normalize"]:
                result /= np.max(np.abs(result))

            if nJob["export"]:
                exportData = np.column_stack((nJob["pathDistances"], result))
                np.savetxt(nJob["label"] + "stage_" + str(nJob["nStages"] - len(self.plotStages)) + ".csv", exportData)

            self.plotter.plotXYData(nJob["pathDistances"], result, nJob["figure"], nJob["axSpec"], nJob)
