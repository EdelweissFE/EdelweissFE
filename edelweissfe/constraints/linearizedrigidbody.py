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

    An implementation of a linearized rigid body constraint for displacements.


    Defining the constraints
    ========================

    2D: 2 constraints per slave node
    1) 0 = ux1 - uxRP + R*dy1
    2) 0 = uy1 - uyRP - R*dx1
    1) 0 = ux2 - uxRP + R*dy2
    2) 0 = uy2 - uyRP - R*dx2
    ...

            | ux1 - uxRP + R*dy1 |   | 0 |
            | uy1 - uyRP - R*dx1 |   | 0 |
    --> g = | ux2 - uxRP + R*dy2 | = | 0 |
            | uy2 - uyRP - R*dx2 |   | 0 |
            | ...                |   | 0 |


    3D: 3 constraints per slave node
    1) 0 = ux1 - uxRP + Rz*dy1 - Ry*dz1
    2) 0 = uy1 - uyRP - Rz*dx1 + Rx*dz1
    3) 0 = uz1 - uzRP - Rx*dy1 + Ry*dx1
    1) 0 = ux2 - uxRP + Rz*dy2 - Ry*dz2
    2) 0 = uy2 - uyRP - Rz*dx2 + Rx*dz2
    3) 0 = uz2 - uzRP - Rx*dy2 + Ry*dx2
    ...

            | ux1 - uxRP + Rz*dy1 - Ry*dz1 |   | 0 |
            | uy1 - uyRP - Rz*dx1 + Rx*dz1 |   | 0 |
            | uz1 - uzRP - Rx*dy1 + Ry*dx1 |   | 0 |
    --> g = | ux2 - uxRP + Rz*dy2 - Ry*dz2 | = | 0 |
            | uy2 - uyRP - Rz*dx2 + Rx*dz2 |   | 0 |
            | uz2 - uzRP - Rx*dy2 + Ry*dx2 |   | 0 |
            | ...                          |   | 0 |



    Create dg/du matrix
    ===================

    2D
            | d_dux1 d_duy1 d_dux2 d_duy2 ... d_duxRP d_duyRP d_dR |
    --------+---------------------------- ... -------------------- +
    dg1_dU  | 1      0      0      0      ... -1       0       dy1 |
    dg2_dU  | 0      1      0      0      ...  0      -1      -dx1 |
    dg1_dU  | 0      0      1      0      ... -1       0       dy2 |
    dg2_dU  | 0      0      0      1      ...  0      -1      -dx2 |
    ...     |......................................................|


    3D
            | d_dux1 d_duy1 d_duz1 d_dux2 d_duy2 d_duz2 ... d_duxRP d_duyRP d_duzRP d_dRx d_dRy d_dRz |
    --------+------------------------------------------ ... ----------------------------------------- +
    dg1_dU  | 1      0      0      0      0      0      ... -1       0       0       0    -dz1   dy1  |
    dg2_dU  | 0      1      0      0      0      0      ...  0      -1       0       dz1   0    -dx1  |
    dg3_dU  | 0      0      1      0      0      0      ...  0       0      -1      -dy1   dx1   0    |
    dg1_dU  | 0      0      0      1      0      0      ... -1       0       0       0    -dz2   dy2  |
    dg2_dU  | 0      0      0      0      1      0      ...  0      -1       0       dz2   0    -dx2  |
    dg3_dU  | 0      0      0      0      0      1      ...  0       0      -1      -dy2   dx2   0    |
    ...     |.........................................................................................|



    Create K matrix
    ===============

    K =     | 0       dg_du.T |
            | dg_du   0       |



    """

    def __init__(self, name, options, model):
        if model.domainSize not in [2, 3]:
            raise WrongDomain("Wrong domain size!")

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

        # slave node set may contain the reference point
        slaveNodeSet = nodeSets[rbNset]

        # RP is removed (if present) and node set is converted to list
        self.slaveNodes = [node for node in slaveNodeSet if not node == self.rp]

        # list of all nodes including RP at end
        self._nodes = self.slaveNodes + [self.rp]

        nSlaves = len(self.slaveNodes)

        self.slaveNodesFields = [["displacement"]] * nSlaves
        self.referencePointFields = [["displacement", "rotation"]]
        self._fieldsOnNodes = self.slaveNodesFields + self.referencePointFields

        nDim = model.domainSize

        nRp = 1

        # one rotational dof for 2D, three for 3D
        nRot = 1 if nDim == 2 else 3

        nConstraints = nSlaves * nDim
        self.nConstraints = nConstraints

        nAffectedDofs = nDim * (nSlaves + nRp) + nRot

        self._nDof = nAffectedDofs + nConstraints

        distances = [s.coordinates - self.rp.coordinates for s in self.slaveNodes]

        dG_dU = np.zeros((nConstraints, nAffectedDofs))

        for i in range(nConstraints):
            # 1) Set the values for d_dux, d_duy and d_duz columns
            # filling diagonal elements with 1
            dG_dU[i, i] = 1

            # 2) Set the values for d_duxRP, d_duyRP and d_duzRP columns
            # This ensures that the block starts at nConstraints and then
            # cycles every nDim rows
            col_index = nConstraints + (i % nDim)

            # make sure we don't go out of bounds
            if col_index < nAffectedDofs - nRot:
                dG_dU[i, col_index] = -1

            # 3) Set the values for the last columns d_dRx, d_dRy, d_dRz
            # index for distances
            index = i // nDim
            # make sure we don't go out of bounds when calling
            # np.array(distances)[index]. It means: if index < number of
            # rows in distances array.
            if index < np.array(distances).shape[0]:
                # 3D
                if nDim == 3:
                    # Accesses the row at the specified index. Unpacks the
                    # x, y, and z components of that row into the variables
                    # dx, dy, and dz, respectively
                    dx, dy, dz = np.array(distances)[index]
                    # Row 0, 3, 6, ...
                    if i % nDim == 0:
                        dG_dU[i, -1] = dy
                        dG_dU[i, -2] = -dz
                    # Row 1, 4, 7, ...
                    elif i % nDim == 1:
                        dG_dU[i, -1] = -dx
                        dG_dU[i, -3] = dz
                    # Row 2, 5, 8, ...
                    elif i % nDim == 2:
                        dG_dU[i, -2] = dx
                        dG_dU[i, -3] = -dy

                # 2D
                elif nDim == 2:
                    dx, dy = np.array(distances)[index]
                    # Row 0, 2, 4, ...
                    if i % nDim == 0:
                        dG_dU[i, -1] = dy
                    # Row 1, 3, 5, ...
                    elif i % nDim == 1:
                        dG_dU[i, -1] = -dx

        self.dG_dU = dG_dU

        K = np.zeros((self._nDof, self._nDof))

        K[0 : dG_dU.shape[1], -dG_dU.shape[0] :] = dG_dU.T
        K[-dG_dU.shape[0] :, 0 : dG_dU.shape[1] :] = dG_dU

        self.K = K

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

        K += self.K
