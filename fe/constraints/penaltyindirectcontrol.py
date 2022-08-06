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
"""
A penalty based constraint used for indirect (displacement) control.
"""

documentation = {
    "field": "the field this constraint acts on",
    "cVector": "the projection vector for the constrained nodes (e.g., CMOD)",
    "constrainedNSet": "the node set for determining the constraint (e.g., CMOD)",
    "loadNSet": "the node set for application of the controlled load",
    "loadVector": "the vector (in correct) dimensions and tensorial order  determining the load",
    "length": "the value of the constraint (e.g., CMOD)",
    "penaltyStiffness": "the stiffness for formulating the constraint",
    "offset": "(optional) a correction value for the computation of the constraint (e.g, initial displacement)",
}

import numpy as np

from fe.config.phenomena import getFieldSize
from fe.utils.misc import stringDict
from fe.utils.exceptions import WrongDomain
from fe.constraints.base.constraintbase import ConstraintBase


class Constraint(ConstraintBase):
    def __init__(self, name, definitionLines, modelInfo):
        super().__init__(name, definitionLines, modelInfo)

        definition = stringDict([e for line in definitionLines for e in line])

        self.theField = definition.get("field", "displacement")

        self.cVector = np.fromstring(definition["cVector"], dtype=np.float, sep=",")
        self.constrainedNSet = modelInfo["nodeSets"][definition["constrainedNSet"]]

        self.loadNSet = modelInfo["nodeSets"][definition["loadNSet"]]
        self.loadVector = np.fromstring(definition["loadVector"], dtype=np.float, sep=",")

        self.penaltyStiffness = float(definition["penaltyStiffness"])
        self.l = np.float(definition["length"])

        if "f(t)" in definition:
            t = sp.symbols("t")
            self.amplitude = sp.lambdify(t, sp.sympify(definition["f(t)"]), "numpy")
        else:
            self.amplitude = lambda x: x

        self.offset = float(definition.get("offset", 0.0))

        self.nodes = self.loadNSet + self.constrainedNSet

        self.fieldsOnNodes = [
            [
                self.theField,
            ]
        ] * len(self.nodes)

        sizeField = getFieldSize(self.theField, modelInfo["domainSize"])

        nDim = modelInfo["domainSize"]

        sizeBlock_loadNodes = nDim * len(self.loadNSet)
        self.startBlock_loadNodes = 0
        self.endBlock_loadNodes = sizeBlock_loadNodes

        sizeBlock_constrainedNodes = nDim * len(self.constrainedNSet)
        self.startBlock_constrainedNodes = self.endBlock_loadNodes
        self.endBlock_constrainedNodes = self.startBlock_constrainedNodes + sizeBlock_constrainedNodes

        self.nDof = self.endBlock_constrainedNodes

        self.active = True

        self.constrainedValue = 0.0

        self.normalizedResidual = np.tile(self.loadVector, len(self.loadNSet))

    def getNumberOfAdditionalNeededScalarVariables(self):
        return 0

    # def assignAdditionalScalarVariables(self, scalarVariables):
    #     pass

    def applyConstraint(self, Un1, dU, PExt, V, increment):

        if self.active == False:
            return

        sBL = self.startBlock_loadNodes
        eBL = self.endBlock_loadNodes
        sBC = self.startBlock_constrainedNodes
        eBC = self.endBlock_constrainedNodes

        U_c = Un1[sBC:eBC]

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment

        L = self.l * self.amplitude(stepProgress)

        cVector = self.cVector

        self.constrainedValue = cVector.dot(U_c)

        loadFactor = self.penaltyStiffness * (self.constrainedValue - self.offset - L)
        dLoadFactor_ddU = self.penaltyStiffness * cVector

        K = V.reshape(self.nDof, self.nDof, order="F")

        t = self.normalizedResidual

        PExt[sBL:eBL] = -t * loadFactor
        K[sBL:eBL, sBC:eBC] = np.outer(t, dLoadFactor_ddU)
