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
# Created on Tue Dec 11 11:21:39 2018

# @author: Matthias Neuner

import numpy as np

from edelweissfe.constraints.base.constraintbase import ConstraintBase
from edelweissfe.utils.misc import convertLinesToStringDictionary

documentation = {
    "nSet": "(slave) node set, which is constrained to the reference point",
    "referencePoint": "(master) reference point",
}


class Constraint(ConstraintBase):
    """
    Geometrically exact rigid body constraint: Constrains a nodeset to a reference point.
    Currently only available for spatialdomain = 3D.
    """

    def __init__(self, name, definitionLines, model):
        super().__init__(name, definitionLines, model)

        # self.name = name
        definition = convertLinesToStringDictionary(definitionLines)

        self.nDim = model.domainSize
        nDim = self.nDim

        if nDim == 2:
            raise Exception("rigid body constraint not yet implemented for 2D")

        rbNset = definition["nSet"]
        nodeSets = model.nodeSets

        if len(nodeSets[definition["referencePoint"]]) > 1:
            raise Exception(
                "node set for reference point '{:}' contains more than one node".format(definition["referencePoint"])
            )

        self.referencePoint = nodeSets[definition["referencePoint"]][0]

        slaveNodeSet = nodeSets[rbNset]  # slave node set may contain the reference point

        # reference point is removed (if present) and node set is converted to list
        self.slaveNodes = [node for node in slaveNodeSet if not node == self.referencePoint]

        nRot = 3
        nSlaves = len(self.slaveNodes)

        self.indicesOfSlaveNodesInP = [[i * nDim + j for j in range(nDim)] for i in range(nSlaves)]
        self.indicesOfRPUinP = [nSlaves * nDim + j for j in range(nDim)]
        self.indicesOfRPPhiInP = [nSlaves * nDim + nDim + j for j in range(nRot)]

        # all nodes

        slaveNodeSet = nodeSets[rbNset]  # slave node set may contain the reference point

        # reference point is removed (if present) and node set is converted to list
        self.slaveNodes = [node for node in slaveNodeSet if not node == self.referencePoint]

        # list of all nodes including RP at end
        self._nodes = self.slaveNodes + [self.referencePoint]

        self.slaveNodesFields = [["displacement"]] * nSlaves
        self.referencePointFields = [["displacement", "rotation"]]
        self._fieldsOnNodes = self.slaveNodesFields + self.referencePointFields

        nConstraints = nSlaves * nDim
        nRp = 1

        self.nDofsOnNodes = nDim * (nSlaves + nRp) + nRot

        self.distancesSlaveNodeRP = [s.coordinates - self.referencePoint.coordinates for s in self.slaveNodes]

        self._nDof = self.nDofsOnNodes + nConstraints
        self.sizeStiffness = self._nDof * self._nDof

        self.nConstraints = nConstraints

        self.nRot = 3

    @property
    def nodes(self) -> list:
        return self._nodes

    @property
    def fieldsOnNodes(self) -> list:
        return self._fieldsOnNodes

    @property
    def nDof(self) -> int:
        return self._nDof

    def update(self, options):
        """No updates are possible for this constraint."""

    def getNumberOfAdditionalNeededScalarVariables(self):
        return self.nConstraints

    def Rz_2D(self, phi, derivative):
        phi = phi + np.pi / 2 * derivative
        return np.array(
            [
                [np.cos(phi), -np.sin(phi)],
                [np.sin(phi), +np.cos(phi)],
            ]
        )

    def Rx_3D(self, phi, derivative):
        phi = phi + np.pi / 2 * derivative
        i = 0.0 if derivative > 0 else 1.0
        return np.array(
            [
                [
                    i,
                    0,
                    0,
                ],
                [0, np.cos(phi), -np.sin(phi)],
                [0, np.sin(phi), +np.cos(phi)],
            ]
        )

    def Ry_3D(self, phi, derivative):
        phi = phi + np.pi / 2 * derivative
        i = 0.0 if derivative > 0 else 1.0
        return np.array(
            [
                [
                    np.cos(phi),
                    0,
                    +np.sin(phi),
                ],
                [0, i, 0],
                [-np.sin(phi), 0, +np.cos(phi)],
            ]
        )

    def Rz_3D(self, phi, derivative):
        phi = phi + np.pi / 2 * derivative
        i = 0.0 if derivative > 0 else 1.0
        return np.array(
            [
                [np.cos(phi), -np.sin(phi), 0],
                [np.sin(phi), +np.cos(phi), 0],
                [
                    0,
                    0,
                    i,
                ],
            ]
        )

    def applyConstraint(self, U_np, dU, PExt, K, timeStep):
        nConstraints = self.nConstraints
        nDim = self.nDim
        nRot = self.nRot

        nU = self._nDof - nConstraints  # nDofs (disp., rot.) without Lagrangian multipliers
        nSlaves = len(self.slaveNodes)

        URp = U_np[self.indicesOfRPUinP]
        PhiRp = U_np[self.indicesOfRPPhiInP]
        Lambdas = U_np[nU:].reshape((nDim, -1), order="F")

        G = np.zeros((nDim, 9))
        H = np.zeros((nDim, 9, 9))

        # dg/dU_Node and # dg/dU_RP
        G[:, 0:nDim] = -np.identity(nDim)
        G[:, nDim : 2 * nDim] = +np.identity(nDim)

        if nDim == 3:
            Rx, Ry, Rz = self.Rx_3D, self.Ry_3D, self.Rz_3D

            RotationMatricesAndDerivatives = [
                [R(phi, derivative) for derivative in range(3)] for R, phi in zip((Rx, Ry, Rz), PhiRp)
            ]
            Rx = RotationMatricesAndDerivatives[0]
            Ry = RotationMatricesAndDerivatives[1]
            Rz = RotationMatricesAndDerivatives[2]

            T = Rz[0] @ Ry[0] @ Rx[0]

            RDerivativeProductsI = (
                Rz[0] @ Ry[0] @ Rx[1],
                Rz[0] @ Ry[1] @ Rx[0],
                Rz[1] @ Ry[0] @ Rx[0],
            )

            RDerivativeProductsII = (
                (
                    Rz[0] @ Ry[0] @ Rx[2],
                    Rz[0] @ Ry[1] @ Rx[1],
                    Rz[1] @ Ry[0] @ Rx[1],
                ),
                (
                    Rz[0] @ Ry[1] @ Rx[1],
                    Rz[0] @ Ry[2] @ Rx[0],
                    Rz[1] @ Ry[1] @ Rx[0],
                ),
                (
                    Rz[1] @ Ry[0] @ Rx[1],
                    Rz[1] @ Ry[1] @ Rx[0],
                    Rz[2] @ Ry[0] @ Rx[0],
                ),
            )

        #        elif nDim == 2:
        #            raise Exception("rigid body constraint not yet implemented for 2D")

        # start and end of Lambda in P
        L0, LF = nU, nU + nDim

        for i in range(nSlaves):
            d0 = self.distancesSlaveNodeRP[i]
            indcsUNode = self.indicesOfSlaveNodesInP[i]
            U_n = U_np[indcsUNode]
            Lambda = Lambdas[:, i]

            g = -d0 - (U_n - URp) + T @ d0

            for i in range(nRot):
                G[:, 2 * nDim + i] = RDerivativeProductsI[i] @ d0
                for j in range(nRot):
                    # only the rot. block is nonzero
                    H[:, 2 * nDim + i, 2 * nDim + j] = RDerivativeProductsII[i][j] @ d0

            indcsU = np.array(indcsUNode + self.indicesOfRPUinP + self.indicesOfRPPhiInP, dtype=int)

            PExt[indcsU] -= Lambda.T @ G
            PExt[L0:LF] -= g

            KUU = K[0:nU, 0:nU]
            KUL = K[0:nU, L0:LF]
            KLU = K[L0:LF, 0:nU]

            # for KUU, only the Phi_RP block is nonzero
            KUU[-nRot:, -nRot:] += np.einsum(
                "i,ijk->jk",
                Lambda,
                H[:, -nRot:, -nRot:],
            )  # L_[i] H_[i,j,k]
            KUL[indcsU, :] += G.T
            KLU[:, indcsU] += G

            L0 += nDim
            LF += nDim
