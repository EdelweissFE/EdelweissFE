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
# Created on Sun May 21 11:34:35 2017

# @author: Matthias Neuner

import numpy as np
import sympy as sp

from edelweissfe.constraints.base.constraintbase import ConstraintBase
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.utils.misc import convertLinesToStringDictionary, strtobool

"""
A penalty based constraint used for indirect (displacement) control.
"""

documentation = {
    "field": "The field this constraint acts on",
    "cVector": "The projection vector for the constrained nodes (e.g., CMOD)",
    "constrainedNSet": "The node set for determining the constraint (e.g., CMOD)",
    "loadNSet": "The node set for application of the controlled load",
    "loadVector": "The vector (in correct) dimensions and tensorial order  determining the load",
    "length": "The value of the constraint (e.g., CMOD)",
    "penaltyStiffness": "The stiffness for formulating the constraint",
    "offset": "(Optional) a correction value for the computation of the constraint (e.g, initial displacement)",
    "normalizeLoad": "(Optional) normalize the applied force per node wrt. the number of nodes, i.e., apply a load irrespective of the total number of nodes in ``loadNSet``",  # noqa: E501
}


class Constraint(ConstraintBase):
    def __init__(self, name, definitionLines, model):
        super().__init__(name, definitionLines, model)

        definition = convertLinesToStringDictionary(definitionLines)

        self.theField = definition.get("field", "displacement")

        self.cVector = np.fromstring(definition["cVector"], dtype=float, sep=",")
        self.constrainedNSet = model.nodeSets[definition["constrainedNSet"]]
        self.loadNSet = model.nodeSets[definition["loadNSet"]]

        self.loadVector = np.fromstring(definition["loadVector"], dtype=float, sep=",")

        # we may normalize in order to end up with an identical load irrespective of the number of nodes
        # in the load node set
        if strtobool(definition.get("normalizeLoad", "True")):
            self.loadVector *= 1.0 / len(self.loadNSet)

        self.penaltyStiffness = float(definition["penaltyStiffness"])
        self.length = float(definition["length"])

        if "f(t)" in definition:
            t = sp.symbols("t")
            self.amplitude = sp.lambdify(t, sp.sympify(definition["f(t)"]), "numpy")
        else:
            self.amplitude = lambda x: x

        self.offset = float(definition.get("offset", 0.0))

        self._nodes = list(self.loadNSet) + list(self.constrainedNSet)

        self._fieldsOnNodes = [
            [
                self.theField,
            ]
        ] * len(self._nodes)

        nDim = model.domainSize

        sizeBlock_loadNodes = nDim * len(self.loadNSet)
        self.startBlock_loadNodes = 0
        self.endBlock_loadNodes = sizeBlock_loadNodes

        sizeBlock_constrainedNodes = nDim * len(self.constrainedNSet)
        self.startBlock_constrainedNodes = self.endBlock_loadNodes
        self.endBlock_constrainedNodes = self.startBlock_constrainedNodes + sizeBlock_constrainedNodes

        self._nDof = self.endBlock_constrainedNodes

        self.active = True

        self.constrainedValue = 0.0

        self.unitResidual = np.tile(self.loadVector, len(self.loadNSet))

    @property
    def nodes(self) -> list:
        return self._nodes

    @property
    def fieldsOnNodes(self) -> list:
        return self._fieldsOnNodes

    @property
    def nDof(self) -> int:
        return self._nDof

    def applyConstraint(self, U_np, dU, PExt, K, timeStep: TimeStep):
        if not self.active:
            return

        sBL = self.startBlock_loadNodes
        eBL = self.endBlock_loadNodes
        sBC = self.startBlock_constrainedNodes
        eBC = self.endBlock_constrainedNodes

        U_c = U_np[sBC:eBC]

        L = self.length * self.amplitude(timeStep.stepProgress)

        cVector = self.cVector

        self.constrainedValue = cVector.dot(U_c)

        loadFactor = self.penaltyStiffness * (self.constrainedValue - self.offset - L)
        dLoadFactor_ddU = self.penaltyStiffness * cVector

        t = self.unitResidual

        PExt[sBL:eBL] = -t * loadFactor
        K[sBL:eBL, sBC:eBC] = np.outer(t, dLoadFactor_ddU)
