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
# Created on Sat Jul 22 14:57:48 2017

# @author: Matthias Neuner
"""
FieldOutputs store all kind of analysis results, 
and are defined via the keyword ``*fieldOutput``.
All fieldOutputs are accessable to all outputmanagers at the end of each 
increment, step and job.
Furthermore, they can be exported to ``*.csv`` files at the end of the analysis job.

ATTENTION: 
    If the results are exported to a .csv file with enabled "saveHistory", 
    the time History is automatically appended to the .csv file"
"""

documentation = {
    "name": "name of the fieldOutput",
    "nSet|elSet|node|element|modelData": "entity, for which the fieldOutput is defined",
    "result": "e.g., U, P, stress, strain ...",
    "quadraturePoint": "for element based fieldOutputs only, integers or slices",
    "f(x)": "(optional), apply math (in each increment)",
    "saveHistory": "(optional), save complete History or only last (increment) result. Default: True (node, element) and False (nSet, elSet)",
    "export": "(optional), export the fieldOutput to a file at the end of the job",
    "f_export(x)": "(optional), apply math on the final result (table)",
}

import numpy as np
from fe.utils.misc import stringDict, strToRange, isInteger
from fe.utils.meshtools import extractNodesFromElementSet
from fe.utils.math import createMathExpression, createModelAccessibleFunction
from fe.utils.elementresultcollector import ElementResultCollector
from fe.utils.dofmanager import DofVector
from fe.journal.journal import Journal
from numpy import ndarray


class FieldOutput:
    """
    Entity of a fieldOutput request.

    Carries the history or the latest result.

    Parameters
    ----------
    modelInfo
        A dictionary containing the model tree.
    definition
        A dictionary containing the definition of the field output.
    journal
        The journal object for logging.
    """

    def __init__(self, modelInfo: dict, definition: dict, journal: Journal):
        self.name = definition["name"]
        self.journal = journal
        self.timeTotal = 0.0

        if "nSet" in definition:
            self.domainType = "nSet"

            self.type = "nodalResult"
            self.nSet = modelInfo["nodeSets"][definition["nSet"]]
            self.nSetName = definition["nSet"]
            self.resultVector = definition["result"]
            self.field = definition["field"]
            self.resultIndices = np.array(
                [[n.fields[self.field]] for n in self.nSet if self.field in n.fields], dtype=np.int
            ).ravel()

        elif "elSet" in definition:
            self.domainType = "elSet"

            if definition["result"] == "U" or definition["result"] == "P":
                self.journal.message(
                    "Converting elSet {:} to a nSet due to requested nodal results".format(definition["elSet"]),
                    self.name,
                )
                # an elset was given, but in fact a nodeset was 'meant': we extract the nodes of the elementset!
                self.type = "nodalResult"
                self.nSet = extractNodesFromElementSet(modelInfo["elementSets"][definition["elSet"]])
                self.nSetName = definition["elSet"]
                self.resultVector = definition["result"]
                self.field = definition["field"]
                self.resultIndices = np.array(
                    [[n.fields[self.field]] for n in self.nSet if self.field in n.fields], dtype=np.int
                ).ravel()
            else:
                # it's really an elSet job
                self.type = "elementResult"
                self.elSet = modelInfo["elementSets"][definition["elSet"]]
                self.elSetName = definition["elSet"]
                self.resultName = definition["result"]

                qp = definition["quadraturePoint"]
                quadraturePoints = strToRange(qp) if not isInteger(qp) else [int(qp)]

                self.quadraturePoints = quadraturePoints

                self.elementResultCollector = ElementResultCollector(self.elSet, quadraturePoints, self.resultName)

        elif "modelData" in definition:
            self.domainType = "model"
            self.type = "modelData"
            self.getCurrentModelDataResult = createModelAccessibleFunction(definition["modelData"], modelInfo=modelInfo)

        else:
            raise Exception("invalid field output requested: " + definition["name"])

        self.appendResults = definition.get("saveHistory", False)
        self.result = [] if self.appendResults else None

        if "f(x)" in definition:
            self.f = createMathExpression(definition["f(x)"])
        else:
            self.f = None

        self.timeHistory = []

        self.export = definition.get("export", False)
        if "f_export(x)" in definition:
            self.f_export = createMathExpression(definition["f_export(x)"])
        else:
            self.f_export = None

    def getLastResult(
        self,
    ) -> ndarray:
        """Get the last result, no matter if the history is stored or not.

        Returns
        -------
        ndarray
            The result array.
        """

        return self.result[-1] if self.appendResults else self.result

    def getResultHistory(
        self,
    ) -> ndarray:
        """Get the history.
        Throws an exception if the history is not stored.

        Returns
        -------
        ndarray
            The result history.
        """

        if not self.appendResults:
            raise Exception(
                "fieldOuput {:} does not save any history; please define it with saveHistory=True!".format(self.name)
            )
        return np.asarray(self.result)

    def getTimeHistory(
        self,
    ) -> ndarray:
        """Get the time history.

        Returns
        -------
        ndarray
            The time history.
        """

        return np.asarray(self.timeHistory)

    def updateResults(self, currentTime: float, U: DofVector, P: DofVector):
        """Update the field output.
        Will use the current solution and reaction vector if result is a nodal result.

        Parameters
        ----------
        currentTime
            The current solution time.
        U
            The current solution vector.
        P
            The current reaction force vector.
        """

        self.timeHistory.append(currentTime)
        incrementResult = None

        if self.type == "nodalResult":
            resVec = U if self.resultVector == "U" else P
            incrementResult = resVec[self.resultIndices].reshape(len(self.nSet), -1)

        elif self.type == "elementResult":
            incrementResult = self.elementResultCollector.getCurrentResults()

        elif self.type == "modelData":
            incrementResult = self.getCurrentModelDataResult()

        if self.f:
            incrementResult = self.f(incrementResult)

        if self.appendResults:
            self.result.append(incrementResult)
        else:
            self.result = incrementResult

    def initializeJob(self, startTime, U, P):
        """Initalize everything. Will also update the results
        based on the proved start time and solution.

        Parameters
        ----------
        startTime
            The initial solution time.
        U
            The initial solution vector.
        P
            The initial reaction force vector.
        """

        self.updateResults(startTime, U, P)

    def initializeStep(self, step, stepActions):
        pass

    def finalizeIncrement(self, U: DofVector, P: DofVector, increment: tuple):
        """Finalize an increment, i.e. store the current results.

        Parameters
        ----------
        U
            The current solution vector.
        P
            The current reaction force vector.
        increment
            The finished time increment.
        """

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        newTime = totalTime + dT
        self.updateResults(newTime, U, P)

    def finalizeStep(
        self,
        U,
        P,
    ):
        pass

    def finalizeJob(
        self,
        U,
        P,
    ):
        """Finalize everything.
        If results are to be exported, it is done now.

        Parameters
        ----------
            U
                The final solution vector.
            P
                The final reaction force vector.
        """

        if self.export:
            res = np.asarray(self.result)
            if self.f_export:
                res = self.f_export(res)
            if res.ndim > 2:
                self.journal.message("Reshaping fieldOutput result for export in .csv file", self.name)
                res = res.reshape((res.shape[0], -1))
            if self.appendResults and res.shape[0] == len(self.timeHistory):
                # we also store the time, if result shape and time history are 'compatible'
                self.journal.message("Adding time history for export in .csv file", self.name)
                time = np.asarray(self.timeHistory).reshape(-1, 1)
                resultTable = np.column_stack((time, res))
            else:
                resultTable = res

            np.savetxt(
                "{:}.csv".format(self.export),
                resultTable,
            )

    def setResults(self, values: ndarray):
        """Modifies a result at it's origin, if possible.
        Throws an exception if not possible.

        Parameters
        ----------
        values
            The values.
        """
        if self.type == "elementResult":
            if self.f:
                raise Exception("cannot set field output for modified results (f(x) != None) !")

            for i, el in enumerate(self.elSet):
                for j, g in enumerate(self.quadraturePoints):
                    theArray = el.getResultArray(self.resultName, g, True)
                    theArray[:] = values[i, j, :]
                    el.acceptLastState()

        else:
            raise Exception("setting field output currently only implemented for element!")

    def __eq__(self, other):
        if type(other) == str:
            return other == self.name
        return self.getLastResult() == other

    def __ne__(self, other):
        return self.getLastResult() != other

    def __lt__(self, other):
        return self.getLastResult() < other

    def __le__(self, other):
        return self.getLastResult() <= other

    def __gt__(self, other):
        return self.getLastResult() > other

    def __ge__(self, other):
        return self.getLastResult() >= other

    def __getitem__(self, index):
        return self.getLastResult()[index]

    def __add__(self, other):
        return self.getLastResult() + other

    def __sub__(self, other):
        return self.getLastResult() - other


class FieldOutputController:
    """
    The central module for managing field outputs, which can be used by output managers.

    Parameters
    ----------
    modelInfo
        A dictionary containing the model tree.
    inputFile
        A dictonary containing the field output definitios.
    journal
        The journal instance for logging.
    """

    def __init__(self, modelInfo: dict, inputFile: dict, journal: Journal):

        self.fieldOutputs = {}

        if not inputFile["*fieldOutput"]:
            return
        definition = inputFile["*fieldOutput"][0]

        for defLine in definition["data"]:
            fpDef = stringDict(defLine)
            self.fieldOutputs[fpDef["name"]] = FieldOutput(modelInfo, fpDef, journal)

    def finalizeIncrement(self, U: DofVector, P: DofVector, increment: tuple):
        """Finalize all field outputs at the end of an increment.

        Parameters
        ----------
        U
            The current solution vector.
        P
            The current raction vector.
        increment
            The time increment.
        """

        for output in self.fieldOutputs.values():
            output.finalizeIncrement(U, P, increment)

    def finalizeStep(
        self,
        U: DofVector,
        P: DofVector,
    ):
        """Finalize all field outputs at the end of a step.

        Parameters
        ----------
        U
            The current solution vector.
        P
            The current reaction vector.
        """

        for output in self.fieldOutputs.values():
            output.finalizeStep(U, P)

    def initializeStep(self, step, stepActions):
        """Initalize an step.

        Parameters
        ----------
        step
            The step information.
        stepActions
            The list of stepActions.
        """

        for output in self.fieldOutputs.values():
            output.initializeStep(step, stepActions)

    def initializeJob(
        self,
        startTime: float,
        U: DofVector,
        P: DofVector,
    ):
        """Initialize all field outputs at the beginning of a job.

        Parameters
        ----------
        startTime
            The initial time.
        U
            The initial solution vector.
        P
            The initial reaction vector.
        """
        for output in self.fieldOutputs.values():
            output.initializeJob(startTime, U, P)

    def finalizeJob(
        self,
        U: DofVector,
        P: DofVector,
    ):
        """Finalize all field outputs at the end of a job.

        Parameters
        ----------
        U
            The final solution vector.
        P
            The final reaction vector.
        """
        for output in self.fieldOutputs.values():
            output.finalizeJob(U, P)
