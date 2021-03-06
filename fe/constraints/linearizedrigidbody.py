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
"""
Created on Sun May 21 11:34:35 2017

@author: Matthias Neuner

Linearized rigid body constraint:
Constrains a nodeset to a reference point.
Currently only available for spatialdomain = 2D.
    
Datalines:
"""


documentation = {
    "nSet": "(slave) node set, which is constrained to the reference point",
    "referencePoint": "(master) reference point",
}

import numpy as np

from fe.utils.misc import stringDict
from fe.utils.exceptions import WrongDomain


class Constraint:
    """linearized Rigid Body constraint"""

    # linear constraints are independent of the solution vector,
    # and, thus, need only be evaluated once (per step)
    linearConstraint = True

    def __init__(self, name, definitionLines, modelInfo):

        if modelInfo["domainSize"] != 2:
            raise WrongDomain("Liniearized Rigid Body is currently only available for 2d domain size")

        self.name = name
        definition = stringDict([e for line in definitionLines for e in line])

        rbNset = definition["nSet"]
        nodeSets = modelInfo["nodeSets"]

        self.rp = nodeSets[definition["referencePoint"]][0]
        self.slaveNodes = nodeSets[rbNset]  # may also contain the RP, doesn't really matter as we remove it

        if self.rp in self.slaveNodes:  # remove the rp from the slave node set
            self.slaveNodes = [s for s in self.slaveNodes if s is not self.rp]

        # all nodes
        self.nodes = self.slaveNodes + [self.rp]

        nSlaves = len(self.slaveNodes)
        self.slaveNodesFields = [["displacement"]] * nSlaves
        self.referencePointFields = [["displacement", "rotation"]]
        self.fieldsOfNodes = self.slaveNodesFields + self.referencePointFields

        nDim = modelInfo["domainSize"]

        nConstraints = nSlaves * 2  # one for distance, one for angle
        nAffectedDofs = nDim * (nSlaves + 1) + 1

        distances = [s.coordinates - self.rp.coordinates for s in self.slaveNodes]
        dMagnitudeSquares = [d @ d for d in distances]

        self.nDof = nAffectedDofs + nConstraints
        self.sizeStiffness = self.nDof * self.nDof

        dG_dU = np.zeros((nConstraints, nAffectedDofs))

        """
               |n1x         n1y         n2x         n2y    ...  rpx         rpy         rpA|
        -------+------------------------------------------ ... ----------------------------+
        dg1_dU |dx1         dy1         0           0      ...  -dx1        -dy1        0  |
        dg2_dU |0           0           dx2         dy2    ...  -dx2        -dy2        0  |
        ...    :...........................................................................:
        dg1A_dU|-d1y/h??     d1x/h??      0           0      ...   d1y/h??     -d1x/h??     -1.|
        dg2A_dU|0           0           -d1y/h??     d1x/h?? ...   d1y/h??      d1x/h??     -1.|
        ...    |...........................................................................|
        """

        # derivatives distance
        for i, d in enumerate(distances):
            dG_dU[i, i * nDim : i * nDim + nDim] = d.T
            dG_dU[i, -(nDim + 1) : -1] = -d.T
        # derivatives angle
        for i, (d, h2) in enumerate(zip(distances, dMagnitudeSquares)):
            x = d.T / h2
            x[1] *= -1
            dG_dU[nSlaves + i, i * nDim : i * nDim + nDim] = x[::-1].T
            dG_dU[nSlaves + i, -(nDim + 1) : -1] = -x[::-1].T
        dG_dU[nSlaves:, -1] = -1

        K = np.zeros((self.nDof, self.nDof))

        """
        K =     |   0       dG_dU.T |
                |   dG_dU   0       |
        """

        K[0 : dG_dU.shape[1], -dG_dU.shape[0] :] = dG_dU.T
        K[-dG_dU.shape[0] :, 0 : dG_dU.shape[1] :] = dG_dU

        self.K = K
        self.dG_dU = dG_dU

        self.additionalGlobalDofIndices = []
        self.nConstraints = nConstraints

    def getNumberOfAdditionalNeededDofs(self):
        return self.nConstraints

    def assignAdditionalGlobalDofIndices(self, additionalGlobalDofIndices):
        self.additionalGlobalDofIndices = additionalGlobalDofIndices

        # TODO: .. move away from constraint, since this should be not handled by a constraint intself
        self.globalDofIndices = np.asarray(
            [
                i
                for node, nodeFields in zip(self.nodes, self.fieldsOfNodes)
                for nodeField in nodeFields
                for i in node.fields[nodeField]
            ]
            + self.additionalGlobalDofIndices
        )

    def applyConstraint(self, Un1, PExt, V, increment):

        LambdaN1 = Un1[-self.nConstraints :]

        PExt[: -self.nConstraints] -= LambdaN1.dot(self.dG_dU)
        #        PExt[- self.nConstraints  : ] -= 0.0

        V += self.K.ravel()
