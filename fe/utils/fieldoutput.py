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
from fe.models.femodel import FEModel
from fe.utils.misc import convertLineToStringDictionary, strToRange, isInteger
from fe.utils.meshtools import extractNodesFromElementSet
from fe.utils.math import createMathExpression, createModelAccessibleFunction
from fe.utils.elementresultcollector import ElementResultCollector
from fe.numerics.dofmanager import DofVector
from fe.models.femodel import FEModel
from fe.journal.journal import Journal
from numpy import ndarray


class _FieldOutput:
    """
    Entity of a fieldOutput request.

    Carries the history or the latest result.

    Parameters
    ----------
    model
        A dictionary containing the model tree.
    definition
        A dictionary containing the definition of the field output.
    journal
        The journal object for logging.
    **kwargs
        The definition for this output.
    """

    def __init__(self, name: str, model: FEModel, journal: Journal, **kwargs: dict):
        self.timeTotal = 0.0
        self.name = name
        self.model = model
        self.journal = journal

        if "nset" in kwargs:
            self.domainType = "nSet"

            self.type = "nodalResult"
            self.nSet = model.nodeSets[kwargs["nset"]]
            self.nSetName = kwargs["nset"]
            self.resultVector = kwargs["result"]
            self.field = kwargs["field"]

        elif "elset" in kwargs:
            self.domainType = "elSet"

            if kwargs["result"] == "U" or kwargs["result"] == "P":
                journal.message(
                    "Converting elSet {:} to a nSet due to requested nodal results".format(kwargs["elset"]),
                    self.name,
                )
                # an elset was given, but in fact a nodeset was 'meant': we extract the nodes of the elementset!
                self.type = "nodalResult"
                self.nSet = extractNodesFromElementSet(model.elementSets[kwargs["elset"]])
                self.nSetName = kwargs["elset"]
                self.resultVector = kwargs["result"]
                self.field = kwargs["field"]
            else:
                # it's really an elSet job
                self.type = "elementResult"
                self.elSet = model.elementSets[kwargs["elset"]]
                self.elSetName = kwargs["elset"]
                self.resultName = kwargs["result"]

                qp = kwargs["quadraturepoint"]
                quadraturePoints = strToRange(qp) if not isInteger(qp) else [int(qp)]

                self.quadraturePoints = quadraturePoints

                self.elementResultCollector = ElementResultCollector(
                    list(self.elSet), quadraturePoints, self.resultName
                )

        elif "modeldata" in kwargs:
            self.domainType = "model"
            self.type = "modelData"
            self.getCurrentModelDataResult = createModelAccessibleFunction(kwargs["modeldata"], model=model)

        else:
            raise Exception("Invalid FieldOutput requested: " + self.name)

        self.appendResults = kwargs.get("savehistory", False)
        self.result = [] if self.appendResults else None

        if "f(x)" in kwargs:
            self.f = createMathExpression(kwargs["f(x)"])
        else:
            self.f = None

        self.timeHistory = []

        self.export = kwargs.get("export", False)
        if "f_export(x)" in kwargs:
            self.f_export = createMathExpression(kwargs["f_export(x)"])
        else:
            self.f_export = None

    def getLastResult(
        self,
    ) -> np.ndarray:
        """Get the last result, no matter if the history is stored or not.

        Returns
        -------
        np.ndarray
            The result array.
        """

        return self.result[-1] if self.appendResults else self.result

    def getResultHistory(
        self,
    ) -> np.ndarray:
        """Get the history.
        Throws an exception if the history is not stored.

        Returns
        -------
        np.ndarray
            The result history.
        """

        if not self.appendResults:
            raise Exception(
                "fieldOuput {:} does not save any history; please define it with saveHistory=True!".format(self.name)
            )
        return np.asarray(self.result)

    def getTimeHistory(
        self,
    ) -> np.ndarray:
        """Get the time history.

        Returns
        -------
        np.ndarray
            The time history.
        """

        return np.asarray(self.timeHistory)

    def updateResults(self, model: FEModel):
        """Update the field output.
        Will use the current solution and reaction vector if result is a nodal result.

        Parameters
        ----------
        model
            The model tree.
        """

        self.timeHistory.append(model.time)
        incrementResult = None

        if self.type == "nodalResult":
            incrementResult = model.nodeFields[self.field].values[self.resultVector]

            if self.nSet and self.nSet != "all":
                incrementResult = np.asarray([incrementResult[n] for n in self.nSet])

        elif self.type == "elementResult":
            incrementResult = self.elementResultCollector.getCurrentResults()

        elif self.type == "modelData":
            incrementResult = self.getCurrentModelDataResult()

        if self.f:
            incrementResult = self.f(incrementResult)

        if self.appendResults:
            self.result.append(incrementResult.copy())
        else:
            self.result = incrementResult

    def initializeJob(
        self,
    ):
        """Initalize everything. Will also update the results
        based on the proved start time and solution.
        """

        self.updateResults(self.model)

    def initializeStep(self, step):
        pass

    def finalizeIncrement(self, increment: tuple):
        """Finalize an increment, i.e. store the current results.

        Parameters
        ----------
        increment
            The finished time increment.
        """

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        self.updateResults(self.model)

    def finalizeStep(
        self,
    ):
        pass

    def finalizeJob(
        self,
    ):
        """Finalize everything.
        If results are to be exported, it is done now.

        """

        if self.export:
            res = np.asarray(self.result)
            if self.f_export:
                res = self.f_export(res)
            if res.ndim > 2:
                journal.message("Reshaping fieldOutput result for export in .csv file", self.name)
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

    def setResults(self, values: np.ndarray):
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
    """

    def __init__(self, model: FEModel, journal: Journal):
        self.model = model
        self.journal = journal
        self.fieldOutputs = {}

    def addFieldOutput(self, name: str, **kwargs: dict):
        """Add a new FieldOutput entry to be computed during the simulation

        Parameters
        ----------
        name
            The name of this FieldOutput.
        model
            The model tree.
        journal
            The Journal instance for logging purposes.
        **kwargs
            The definition of the FieldOutput
        """
        if name in self.fieldOutputs:
            raise Exception("FieldOutput {:} already exists!".format(name))
        self.fieldOutputs[name] = _FieldOutput(name, self.model, self.journal, **kwargs)

    def initializeJob(self):
        for fieldOutput in self.fieldOutputs.values():
            fieldOutput.initializeJob()

    def finalizeIncrement(self, increment: tuple):
        """Finalize all field outputs at the end of an increment.

        Parameters
        ----------
        model
            The model tree.
        increment
            The time increment.
        """

        for output in self.fieldOutputs.values():
            output.finalizeIncrement(increment)

    def finalizeStep(
        self,
    ):
        """Finalize all field outputs at the end of a step."""

        for output in self.fieldOutputs.values():
            output.finalizeStep()

    def initializeStep(self, step):
        """Initalize an step.

        Parameters
        ----------
        step
            The step information.
        """

        for output in self.fieldOutputs.values():
            output.initializeStep(step)

    def finalizeJob(
        self,
    ):
        """Finalize all field outputs at the end of a job."""
        for output in self.fieldOutputs.values():
            output.finalizeJob()
