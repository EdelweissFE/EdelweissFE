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

from edelweissfe.constraints.base.constraintbase import ConstraintBase
from edelweissfe.utils.exceptions import WrongDomain
from edelweissfe.utils.misc import convertLinesToStringDictionary

documentation = {
    "nSet": "(Slave) node set, which is constrained to the reference point",
    "referencePoint": "(Master) reference point",
}


class Constraint(ConstraintBase):
    """Linearized Rigid Body constraint.

    A very simple implementation of a linearized rigid body constraint for displacements, currently only in 2D.
    """

    def __init__(self, name, options, model):
        if model.domainSize != 2:
            raise WrongDomain("Linearized rigid body constraint is currently only available for 2d domain size")

        # self.name = name
        definition = convertLinesToStringDictionary(options)

        rbNset = definition["nSet"]
        nodeSets = model.nodeSets

        rpNodeSet = nodeSets[definition["referencePoint"]]

        if len(rpNodeSet) > 1:
            raise Exception(
                "node set for reference point '{:}' contains more than one node".format(definition["referencePoint"])
            )

        self.rp = rpNodeSet[0]

        slaveNodeSet = nodeSets[rbNset]  # slave node set may contain the reference point

        # reference point is removed (if present) and node set is converted to list
        self.slaveNodes = [node for node in slaveNodeSet if not node == self.rp]

        # list of all nodes including RP at end
        self._nodes = self.slaveNodes + [self.rp]

        nSlaves = len(self.slaveNodes)
        self.slaveNodesFields = [["displacement"]] * nSlaves
        self.referencePointFields = [["displacement", "rotation"]]
        self._fieldsOnNodes = self.slaveNodesFields + self.referencePointFields

        nDim = model.domainSize

        nConstraints = nSlaves * 2  # one for distance, one for angle
        nAffectedDofs = nDim * (nSlaves + 1) + 1

        distances = [s.coordinates - self.rp.coordinates for s in self.slaveNodes]
        dMagnitudeSquares = [d @ d for d in distances]

        self._nDof = nAffectedDofs + nConstraints
        self.sizeStiffness = self._nDof * self._nDof

        dG_dU = np.zeros((nConstraints, nAffectedDofs))

        """
               |n1x         n1y         n2x         n2y    ...  rpx         rpy         rpA|
        -------+------------------------------------------ ... ----------------------------+
        dg1_dU |dx1         dy1         0           0      ...  -dx1        -dy1        0  |
        dg2_dU |0           0           dx2         dy2    ...  -dx2        -dy2        0  |
        ...    :...........................................................................:
        dg1A_dU|-d1y/h²     d1x/h²      0           0      ...   d1y/h²     -d1x/h²     -1.|
        dg2A_dU|0           0           -d1y/h²     d1x/h² ...   d1y/h²      d1x/h²     -1.|
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

        K = np.zeros((self._nDof, self._nDof))

        """
        K =     |   0       dG_dU.T |
                |   dG_dU   0       |
        """

        K[0 : dG_dU.shape[1], -dG_dU.shape[0] :] = dG_dU.T
        K[-dG_dU.shape[0] :, 0 : dG_dU.shape[1] :] = dG_dU

        self.K = K
        self.dG_dU = dG_dU

        self.nConstraints = nConstraints

    @property
    def nodes(self) -> list:
        return self._nodes

    @property
    def fieldsOnNodes(self) -> list:
        return self._fieldsOnNodes

    @property
    def nDof(self) -> int:
        return self._nDof

    def getNumberOfAdditionalNeededScalarVariables(self):
        return self.nConstraints

    def applyConstraint(self, U_np, dU, PExt, K, timeStep):
        LambdaN1 = U_np[-self.nConstraints :]

        PExt[: -self.nConstraints] -= LambdaN1.dot(self.dG_dU)
        #        PExt[- self.nConstraints  : ] -= 0.0

        K += self.K
