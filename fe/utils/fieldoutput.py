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

## @author: Matthias Neuner
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
from fe.fields.nodefield import NodeField, NodeFieldSubset
from fe.sets.elementset import ElementSet


class _FieldOutputBase:
    """
    Entity of a fieldOutput request.

    Carries the history or the latest result.

    Parameters
    ----------
    name
        The name of this FieldOutput.
    model
        A dictionary containing the model tree.
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
        self.appendResults = kwargs.get("saveHistory", False)
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

    def _applyResultsPipleline(self, result):
        """Apply the pipeline of operations onto the results.
        Called by inheriting classes.

        Parameters
        ----------
        model
            The model tree.
        """

        self.timeHistory.append(self.model.time)

        if self.f:
            result = self.f(result)

        if self.appendResults:
            self.result.append(result.copy())
        else:
            self.result = result

    def writeLastResult(self):
        """Update file output.

        Parameters
        ----------
        model
            The model tree.
        """
        res = np.asarray(self.result)
        if self.f_export:
            res = self.f_export(res)

        if res.ndim > 2:
            self.journal.message("Reshaping fieldOutput result for export in .csv file", self.name)
            res = res.reshape((res.shape[0], -1))

        with open(f"{self.export}.csv", "a") as f:
            np.savetxt(
                f,
                np.hstack(
                    (
                        self.timeHistory[-1],
                        res[-1,],
                    )
                ).reshape((1, -1)),
            )

    def initializeJob(self):
        """Initalize everything. Will also update the results
        based on the proved start time and solution.
        """

        self.updateResults(self.model)

        if self.export:
            f = open(f"{self.export}.csv", "w")
            f.close()

    def initializeStep(self, step):
        if self.export:
            self.writeLastResult()

    def finalizeIncrement(
        self,
    ):
        """Finalize an increment, i.e. store the current results."""
        self.updateResults(self.model)

        if self.export:
            self.writeLastResult()

    def finalizeStep(
        self,
    ):
        pass

    def finalizeJob(
        self,
    ):
        pass

    def setResults(self, values: np.ndarray):
        """Modifies a result at it's origin, if possible.
        Throws an exception if not possible.

        Parameters
        ----------
        values
            The values.
        """
        raise Exception("setting field output currently not implemented for this type of output!")

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


class NodeFieldOutput(_FieldOutputBase):
    """
    This is a Node based FieldOutput.
    It operates on NodeFields.

    Parameters
    ----------
    name
        The name of this FieldOutput.
    nodeField
        The NodeField, on which this FieldOutput operates.
    result
        The name of the result entry in the NodeField.
    model
        The model tree instance.
    journal
        The journal object for logging.
    **kwargs
        The definition for this output.
    """

    def __init__(self, name: str, nodeField, result: str, model: FEModel, journal: Journal, **kwargs: dict):
        self.entry = result

        self._nodeField = nodeField

        # if isinstance(nodeField, NodeFieldSubset):
        self.associatedSet = nodeField.associatedSet

        super().__init__(name, model, journal, **kwargs)

    def updateResults(self, model: FEModel):
        """Update the field output.
        Will use the current solution and reaction vector if result is a nodal result.

        Parameters
        ----------
        model
            The model tree.
        """

        result = self._nodeField[self.entry]

        return super()._applyResultsPipleline(result)


class ElementFieldOutput(_FieldOutputBase):
    """
    This is a Element based FieldOutput.
    It operates on ElementSets.

    Parameters
    ----------
    name
        The name of this FieldOutput.
    elSet
        The ElementSet on which this FieldOutput operates.
    resultName
        The name of the result entry in the :class:`ElementBase.
    model
        The model tree instance.
    journal
        The journal object for logging.
    **kwargs
        The definition for this output.
    """

    def __init__(self, name: str, elSet, resultName: str, model: FEModel, journal: Journal, **kwargs: dict):
        self.associatedSet = elSet

        self.resultName = resultName
        qp = kwargs["quadraturePoint"]
        self.quadraturePoints = strToRange(qp) if not isInteger(qp) else [int(qp)]

        self.elementResultCollector = ElementResultCollector(
            list(self.associatedSet), self.quadraturePoints, self.resultName
        )

        super().__init__(name, model, journal, **kwargs)

    def updateResults(self, model: FEModel):
        """Update the field output.
        Will use the current solution and reaction vector if result is a nodal result.

        Parameters
        ----------
        model
            The model tree.
        """

        result = self.elementResultCollector.getCurrentResults()

        super()._applyResultsPipleline(result)

    def setResults(self, values: np.ndarray):
        """Modifies a result at it's origin, if possible.
        Throws an exception if not possible.

        Parameters
        ----------
        values
            The values.
        """
        if self.f:
            raise Exception("cannot set field output for modified results (f(x) != None) !")

        for i, el in enumerate(self.associatedSet):
            for j, g in enumerate(self.quadraturePoints):
                theArray = el.getResultArray(self.resultName, g, True)
                theArray[:] = values[i, j, :]
                el.acceptLastState()


class ExpressionFieldOutput(_FieldOutputBase):
    """
    This is a Node based FieldOutput.
    It operates on NodeFields.

    Parameters
    ----------
    name
        The name of this FieldOutput.
    nodeField
        The NodeField, on which this FieldOutput operates.
    result
        The name of the result entry in the NodeField.
    model
        The model tree instance.
    journal
        The journal object for logging.
    **kwargs
        The definition for this output.
    """

    def __init__(self, name: str, model: FEModel, journal: Journal, **kwargs: dict):
        if "nSet" in kwargs:
            self.associatedSet = kwargs.pop("nSet")
        elif "elSet" in kwargs:
            self.associatedSet = kwargs.pop("elSet")
        else:
            raise Exception("All FieldOuputs must be associated with a set!")

        self.theExpression = createModelAccessibleFunction(kwargs["expression"], model)

        super().__init__(name, model, journal, **kwargs)

    def updateResults(self, model: FEModel):
        """Update the field output.
        Will use the current solution and reaction vector if result is a nodal result.

        Parameters
        ----------
        model
            The model tree.
        """

        result = np.reshape(np.asarray(self.theExpression()), (len(self.associatedSet), -1))

        return super()._applyResultsPipleline(result)


class FieldOutputController:
    """
    The central module for managing field outputs, which can be used by output managers.
    """

    def __init__(self, model: FEModel, journal: Journal):
        self.model = model
        self.journal = journal
        self.fieldOutputs = {}

    def addExpressionFieldOutput(self, name: str, **kwargs: dict):
        """Add a new FieldOutput entry to be computed during the simulation

        Parameters
        ----------
        name
            The name of this FieldOutput.
        nodeField
            The :class:`NodeField, on which this FieldOutput should operate.
        resultName
            The name of the result entry in the :class:`NodeField
        **kwargs
            Further definitions of the FieldOutput
        """

        if name in self.fieldOutputs:
            raise Exception("FieldOutput {:} already exists!".format(name))

        self.fieldOutputs[name] = ExpressionFieldOutput(name, self.model, self.journal, **kwargs)

    def addPerNodeFieldOutput(self, name: str, nodeField: NodeField, result: str, **kwargs: dict):
        """Add a new FieldOutput entry to be computed during the simulation

        Parameters
        ----------
        name
            The name of this FieldOutput.
        nodeField
            The :class:`NodeField, on which this FieldOutput should operate.
        resultName
            The name of the result entry in the :class:`NodeField
        **kwargs
            Further definitions of the FieldOutput
        """

        if name in self.fieldOutputs:
            raise Exception("FieldOutput {:} already exists!".format(name))

        self.fieldOutputs[name] = NodeFieldOutput(name, nodeField, result, self.model, self.journal, **kwargs)

    def addPerElementFieldOutput(self, name: str, elSet: ElementSet, resultName: str, **kwargs: dict):
        """Add a new FieldOutput entry to be computed during the simulation

        Parameters
        ----------
        name
            The name of this FieldOutput.
        elSet
            The :class:`ElementSet on which this FieldOutput should operate.
        resultname
            The name of the result, which is provided by the Elements in the :class:`ElementSet.
        journal
            The :class:`Journal instance for logging purposes.
        **kwargs
            Further definition of the FieldOutput
        """

        if name in self.fieldOutputs:
            raise Exception("FieldOutput {:} already exists!".format(name))

        self.fieldOutputs[name] = ElementFieldOutput(name, elSet, resultName, self.model, self.journal, **kwargs)

    def initializeJob(self):
        for fieldOutput in self.fieldOutputs.values():
            fieldOutput.initializeJob()

    def finalizeIncrement(
        self,
    ):
        """Finalize all field outputs at the end of an increment."""

        for output in self.fieldOutputs.values():
            output.finalizeIncrement()

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
