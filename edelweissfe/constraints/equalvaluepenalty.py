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
# Created on Fri Sep 9 11:34:35 2022

# @author: Matthias Neuner

import numpy as np

from edelweissfe.config.phenomena import getFieldSize
from edelweissfe.constraints.base.constraintbase import ConstraintBase
from edelweissfe.models.femodel import FEModel
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.utils.misc import convertLinesToStringDictionary

"""
A penalty based constraint used for constraining nodal values
of a node set to be equal.
"""

documentation = {
    "field": "The field this constraint acts on.",
    "component": "The component of the field.",
    "penalty": "The numerical penalty value.",
    "nSet": "The node set to be constrained.",
}


class Constraint(ConstraintBase):
    def __init__(self, name: str, definitionLines: list, model: FEModel):
        super().__init__(name, definitionLines, model)

        definition = convertLinesToStringDictionary(definitionLines)

        theField = definition["field"]
        self.sizeField = getFieldSize(theField, model.domainSize)
        self.component = int(definition["component"])
        self.penalty = float(definition["penalty"])
        self._nodes = model.nodeSets[definition["nSet"]]
        self._nNodes = len(self._nodes)
        self._nDof = self.sizeField * self._nNodes

        self.indices_component = slice(self.component, self._nDof + self.component, self.sizeField)

        self._fieldsOnNodes = [
            [
                theField,
            ]
        ] * self._nNodes

        self.active = True

    @property
    def nodes(self) -> list:
        return self._nodes

    @property
    def fieldsOnNodes(self) -> list:
        return self._fieldsOnNodes

    @property
    def nDof(self) -> int:
        return self._nDof

    def applyConstraint(
        self,
        U_np: np.ndarray,
        dU: np.ndarray,
        PExt: np.ndarray,
        K: np.ndarray,
        timeStep: TimeStep,
    ):
        if not self.active:
            return

        values = U_np[self.indices_component]

        mean = np.sum(values) / self._nNodes

        PExt[self.indices_component] -= self.penalty * (values - mean)

        diag = np.diag(K)
        diag.setflags(write=True)  # bug in numpy
        diag[self.indices_component] += self.penalty
        K[self.indices_component, self.indices_component] += -self.penalty * 1.0 / self._nNodes
