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
#  Daniel Reitmair daniel.reitmair@uibk.ac.at
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


def computeJacobian(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, x: np.ndarray, nInt: int, nNodes: int, dim: int):
    """Get the Jacobi matrix for the element calculation.

    Parameters
    ----------
    xi
        Local coordinate xi.
    eta
        Local coordinate eta.
    z
        Local coordinate zeta.
    x
        Global coordinates of the element points.
    nInt
        Number of integration points.
    nNodes
        Number of nodes the element has.
    dim
        Dimension the element has."""

    if dim == 2:
        if nNodes == 4:
            return _J2D4(xi, eta, x, nInt)
        elif nNodes == 8:
            return _J2D8(xi, eta, x, nInt)
    elif dim == 3:
        if nNodes == 8:
            return _J3D8(xi, eta, z, x, nInt)
        elif nNodes == 20:
            return _J3D20(xi, eta, z, x, nInt)


def computeBOperator(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, x: np.ndarray, nInt: int, nNodes: int, dim: int):
    """Get the B operator for the element calculation.

    Parameters
    ----------
    xi
        Local coordinate xi.
    eta
        Local coordinate eta.
    z
        Local coordinate zeta.
    x
        Global coordinates of the element points.
    nInt
        Number of integration points.
    nNodes
        Number of nodes the element has.
    dim
        Dimension the element has."""

    if dim == 2:
        if nNodes == 4:
            return _B2D4(xi, eta, x, nInt)
        elif nNodes == 8:
            return _B2D8(xi, eta, x, nInt)
    elif dim == 3:
        if nNodes == 8:
            return _B3D8(xi, eta, z, x, nInt)
        elif nNodes == 20:
            return _B3D20(xi, eta, z, x, nInt)


def computeNOperator(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, nInt: int, nNodes: int, dim: int):
    """Get the N operator containing the shape functions.

    Parameters
    ----------
    xi
        Local coordinates xi for the integration points.
    eta
        Local coordinates eta for the integration points.
    z
        Local coordinates zeta for the integration points.
    nInt
        Number of integration points.
    nNodes
        Number of nodes the element has.
    dim
        Dimension the element has.
    """

    N = np.zeros([nInt, nNodes])
    if dim == 2:
        if nNodes == 4:  # Quad4
            for i in range(nInt):
                N[i] = (
                    1
                    / 4
                    * np.array(
                        [
                            (1 - xi[i]) * (1 - eta[i]),
                            (1 + xi[i]) * (1 - eta[i]),
                            (1 + xi[i]) * (1 + eta[i]),
                            (1 - xi[i]) * (1 + eta[i]),
                        ]
                    )
                )
        elif nNodes == 8:  # Quad8
            for i in range(nInt):
                N[i] = (
                    1
                    / 4
                    * np.array(
                        [
                            (1 - xi[i]) * (1 - eta[i]) * (-xi[i] - eta[i] - 1),
                            (1 + xi[i]) * (1 - eta[i]) * (xi[i] - eta[i] - 1),
                            (1 + xi[i]) * (1 + eta[i]) * (xi[i] + eta[i] - 1),
                            (1 - xi[i]) * (1 + eta[i]) * (-xi[i] + eta[i] - 1),
                            2 * (1 - xi[i] ** 2) * (1 - eta[i]),
                            2 * (1 + xi[i]) * (1 - eta[i] ** 2),
                            2 * (1 - xi[i] ** 2) * (1 + eta[i]),
                            2 * (1 - xi[i]) * (1 - eta[i] ** 2),
                        ]
                    )
                )
    elif dim == 3:
        if nNodes == 8:  # Hexa8
            for i in range(nInt):
                N[i] = (
                    1
                    / 8
                    * np.array(
                        [
                            (1 - xi[i]) * (1 - eta[i]) * (1 - z[i]),
                            (1 - xi[i]) * (1 - eta[i]) * (1 + z[i]),
                            (1 + xi[i]) * (1 - eta[i]) * (1 + z[i]),
                            (1 + xi[i]) * (1 - eta[i]) * (1 - z[i]),
                            (1 - xi[i]) * (1 + eta[i]) * (1 - z[i]),
                            (1 - xi[i]) * (1 + eta[i]) * (1 + z[i]),
                            (1 + xi[i]) * (1 + eta[i]) * (1 + z[i]),
                            (1 + xi[i]) * (1 + eta[i]) * (1 - z[i]),
                        ]
                    )
                )
        elif nNodes == 20:  # Hexa20
            for i in range(nInt):
                N[i] = (
                    1
                    / 8
                    * np.array(
                        [
                            (1 - xi[i]) * (1 - eta[i]) * (1 - z[i]) * (-xi[i] - eta[i] - z[i] - 2),
                            (1 - xi[i]) * (1 - eta[i]) * (1 + z[i]) * (-xi[i] - eta[i] + z[i] - 2),
                            (1 + xi[i]) * (1 - eta[i]) * (1 + z[i]) * (xi[i] - eta[i] + z[i] - 2),
                            (1 + xi[i]) * (1 - eta[i]) * (1 - z[i]) * (xi[i] - eta[i] - z[i] - 2),
                            (1 - xi[i]) * (1 + eta[i]) * (1 - z[i]) * (-xi[i] + eta[i] - z[i] - 2),
                            (1 - xi[i]) * (1 + eta[i]) * (1 + z[i]) * (-xi[i] + eta[i] + z[i] - 2),
                            (1 + xi[i]) * (1 + eta[i]) * (1 + z[i]) * (xi[i] + eta[i] + z[i] - 2),
                            (1 + xi[i]) * (1 + eta[i]) * (1 - z[i]) * (xi[i] + eta[i] - z[i] - 2),
                            2 * (1 - xi[i]) * (1 - eta[i]) * (1 - z[i] ** 2),
                            2 * (1 - xi[i] ** 2) * (1 - eta[i]) * (1 + z[i]),
                            2 * (1 + xi[i]) * (1 - eta[i]) * (1 - z[i] ** 2),
                            2 * (1 - xi[i] ** 2) * (1 - eta[i]) * (1 - z[i]),
                            2 * (1 - xi[i]) * (1 + eta[i]) * (1 - z[i] ** 2),
                            2 * (1 - xi[i] ** 2) * (1 + eta[i]) * (1 + z[i]),
                            2 * (1 + xi[i]) * (1 + eta[i]) * (1 - z[i] ** 2),
                            2 * (1 - xi[i] ** 2) * (1 + eta[i]) * (1 - z[i]),
                            2 * (1 - xi[i]) * (1 - eta[i] ** 2) * (1 - z[i]),
                            2 * (1 - xi[i]) * (1 - eta[i] ** 2) * (1 + z[i]),
                            2 * (1 + xi[i]) * (1 - eta[i] ** 2) * (1 + z[i]),
                            2 * (1 + xi[i]) * (1 - eta[i] ** 2) * (1 - z[i]),
                        ]
                    )
                )
    return N


def _J2D4(xi: np.ndarray, eta: np.ndarray, x: np.ndarray, nInt: int):
    """Get the Jacobi matrix for a Quad4 element."""

    J = np.zeros([nInt, 2, 2])
    # calc all parameters for the X and Y functions (Q4)
    A = np.array([[1, -1, -1, 1], [1, 1, -1, -1], [1, 1, 1, 1], [1, -1, 1, -1]])
    invA = np.linalg.inv(A)
    # calculate parameters
    ax = invA @ np.transpose(x[0])
    ay = invA @ np.transpose(x[1])
    for i in range(0, nInt):  # for all Gauss points (N in total)
        # [J] Jacobi matrix (only Q4)
        J[i] = np.array(
            [
                [ax[1] + ax[3] * xi[i], ay[1] + ay[3] * xi[i]],
                [ax[2] + ax[3] * eta[i], ay[2] + ay[3] * eta[i]],
            ]
        )
    return J


def _J2D8(xi: np.ndarray, eta: np.ndarray, x: np.ndarray, nInt: int):
    """Get the Jacobi matrix for a Quad8 element."""

    J = np.zeros([nInt, 2, 2])
    # calc all parameters for the X and Y functions (Q8)
    A = np.array(
        [
            [1, -1, -1, 1, 1, 1, -1, -1],
            [1, 1, -1, -1, 1, 1, -1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, -1, 1, -1, 1, 1, 1, -1],
            [1, 0, -1, 0, 0, 1, 0, 0],
            [1, 1, 0, 0, 1, 0, 0, 0],
            [1, 0, 1, 0, 0, 1, 0, 0],
            [1, -1, 0, 0, 1, 0, 0, 0],
        ]
    )
    invA = np.linalg.inv(A)
    # calculate parameters
    ax = invA @ np.transpose(x[0])
    ay = invA @ np.transpose(x[1])
    for i in range(0, nInt):  # for all Gauss points (N in total)
        # [J] Jacobi matrix for Q8
        J[i] = np.array(
            [
                [
                    ax[1] + ax[3] * xi[i] + 2 * ax[4] * eta[i] + 2 * ax[6] * eta[i] * xi[i] + ax[7] * xi[i] ** 2,
                    ay[1] + ay[3] * xi[i] + 2 * ay[4] * eta[i] + 2 * ay[6] * eta[i] * xi[i] + ay[7] * xi[i] ** 2,
                ],
                [
                    ax[2] + ax[3] * eta[i] + 2 * ax[5] * xi[i] + ax[6] * eta[i] ** 2 + 2 * ax[7] * eta[i] * xi[i],
                    ay[2] + ay[3] * eta[i] + 2 * ay[5] * xi[i] + ay[6] * eta[i] ** 2 + 2 * ay[7] * eta[i] * xi[i],
                ],
            ]
        )
    return J


def _J3D8(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, x: np.ndarray, nInt: int):
    """Get the Jacobi matrix for a Hexa8 element."""

    J = np.zeros([nInt, 3, 3])
    # calc all parameters for the X and Y functions (H8)
    A = np.array(
        [
            [1, -1, -1, -1, 1, 1, 1, -1],
            [1, -1, -1, 1, 1, -1, -1, 1],
            [1, 1, -1, 1, -1, -1, 1, -1],
            [1, 1, -1, -1, -1, 1, -1, 1],
            [1, -1, 1, -1, -1, -1, 1, 1],
            [1, -1, 1, 1, -1, 1, -1, -1],
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, -1, 1, -1, -1, -1],
        ]
    )
    invA = np.linalg.inv(A)
    ax = invA @ np.transpose(x[0])
    ay = invA @ np.transpose(x[1])
    az = invA @ np.transpose(x[2])
    for i in range(0, nInt):  # for all Gauss points (N in total)
        # [J] Jacobi matrix (only H8)
        J[i] = np.array(
            [
                [
                    ax[1] + ax[4] * xi[i] + ax[6] * z[i] + ax[7] * xi[i] * z[i],
                    ay[1] + ay[4] * xi[i] + ay[6] * z[i] + ay[7] * xi[i] * z[i],
                    az[1] + az[4] * xi[i] + az[6] * z[i] + az[7] * xi[i] * z[i],
                ],
                [
                    ax[2] + ax[4] * eta[i] + ax[5] * z[i] + ax[7] * eta[i] * z[i],
                    ay[2] + ay[4] * eta[i] + ay[5] * z[i] + ay[7] * eta[i] * z[i],
                    az[2] + az[4] * eta[i] + az[5] * z[i] + az[7] * eta[i] * z[i],
                ],
                [
                    ax[3] + ax[5] * xi[i] + ax[6] * eta[i] + ax[7] * eta[i] * xi[i],
                    ay[3] + ay[5] * xi[i] + ay[6] * eta[i] + ay[7] * eta[i] * xi[i],
                    az[3] + az[5] * xi[i] + az[6] * eta[i] + az[7] * eta[i] * xi[i],
                ],
            ]
        )
    return J


def _J3D20(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, x: np.ndarray, nInt: int):
    """Get the Jacobi matrix for a Hexa20 element."""

    J = np.zeros([nInt, 3, 3])
    # calc all parameters for the X and Y functions (H20)
    A = np.array(
        [
            [1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1],
            [1, -1, -1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1],
            [1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1],
            [1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1],
            [1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1],
            [1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, -1, 1, -1, -1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, -1, 1, -1, -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1],
            [1, -1, -1, 0, 1, 0, 0, 1, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, -1, 1, 0, -1, 0, 0, 1, 1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
            [1, 1, -1, 0, -1, 0, 0, 1, 1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, -1, -1, 0, 1, 0, 0, 1, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0],
            [1, -1, 1, 0, -1, 0, 0, 1, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0],
            [1, -1, 0, -1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0],
            [1, -1, 0, 1, 0, 0, -1, 1, 0, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0],
            [1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
            [1, 1, 0, -1, 0, 0, -1, 1, 0, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0],
        ]
    )
    invA = np.linalg.inv(A)
    ax = invA @ np.transpose(x[0])
    ay = invA @ np.transpose(x[1])
    az = invA @ np.transpose(x[2])
    for i in range(0, nInt):  # for all Gauss points (N in total)
        # [J] Jacobi matrix (only H8)
        J[i] = np.array(
            [
                [
                    ax[1]
                    + ax[4] * xi[i]
                    + ax[6] * z[i]
                    + 2 * ax[7] * eta[i]
                    + 2 * ax[10] * eta[i] * xi[i]
                    + ax[11] * xi[i] ** 2
                    + ax[14] * z[i] ** 2
                    + 2 * ax[15] * z[i] * eta[i]
                    + ax[16] * xi[i] * z[i]
                    + 2 * ax[17] * eta[i] * xi[i] * z[i]
                    + ax[18] * xi[i] ** 2 * z[i]
                    + ax[19] * xi[i] * z[i] ** 2,
                    ay[1]
                    + ay[4] * xi[i]
                    + ay[6] * z[i]
                    + 2 * ay[7] * eta[i]
                    + 2 * ay[10] * eta[i] * xi[i]
                    + ay[11] * xi[i] ** 2
                    + ay[14] * z[i] ** 2
                    + 2 * ay[15] * z[i] * eta[i]
                    + ay[16] * xi[i] * z[i]
                    + 2 * ay[17] * eta[i] * xi[i] * z[i]
                    + ay[18] * xi[i] ** 2 * z[i]
                    + ay[19] * xi[i] * z[i] ** 2,
                    az[1]
                    + az[4] * xi[i]
                    + az[6] * z[i]
                    + 2 * az[7] * eta[i]
                    + 2 * az[10] * eta[i] * xi[i]
                    + az[11] * xi[i] ** 2
                    + az[14] * z[i] ** 2
                    + 2 * az[15] * z[i] * eta[i]
                    + az[16] * xi[i] * z[i]
                    + 2 * az[17] * eta[i] * xi[i] * z[i]
                    + az[18] * xi[i] ** 2 * z[i]
                    + az[19] * xi[i] * z[i] ** 2,
                ],
                [
                    ax[2]
                    + ax[4] * eta[i]
                    + ax[5] * z[i]
                    + 2 * ax[8] * xi[i]
                    + ax[10] * eta[i] ** 2
                    + 2 * ax[11] * eta[i] * xi[i]
                    + 2 * ax[12] * xi[i] * z[i]
                    + ax[13] * z[i] ** 2
                    + ax[16] * eta[i] * z[i]
                    + ax[17] * eta[i] ** 2 * z[i]
                    + 2 * ax[18] * eta[i] * xi[i] * z[i]
                    + ax[19] * eta[i] * z[i] ** 2,
                    ay[2]
                    + ay[4] * eta[i]
                    + ay[5] * z[i]
                    + 2 * ay[8] * xi[i]
                    + ay[10] * eta[i] ** 2
                    + 2 * ay[11] * eta[i] * xi[i]
                    + 2 * ay[12] * xi[i] * z[i]
                    + ay[13] * z[i] ** 2
                    + ay[16] * eta[i] * z[i]
                    + ay[17] * eta[i] ** 2 * z[i]
                    + 2 * ay[18] * eta[i] * xi[i] * z[i]
                    + ay[19] * eta[i] * z[i] ** 2,
                    az[2]
                    + az[4] * eta[i]
                    + az[5] * z[i]
                    + 2 * az[8] * xi[i]
                    + az[10] * eta[i] ** 2
                    + 2 * az[11] * eta[i] * xi[i]
                    + 2 * az[12] * xi[i] * z[i]
                    + az[13] * z[i] ** 2
                    + az[16] * eta[i] * z[i]
                    + az[17] * eta[i] ** 2 * z[i]
                    + 2 * az[18] * eta[i] * xi[i] * z[i]
                    + az[19] * eta[i] * z[i] ** 2,
                ],
                [
                    ax[3]
                    + ax[5] * z[i]
                    + ax[6] * eta[i]
                    + 2 * ax[9] * z[i]
                    + ax[12] * xi[i] ** 2
                    + 2 * ax[13] * xi[i] * z[i]
                    + 2 * ax[14] * z[i] * eta[i]
                    + ax[15] * eta[i] ** 2
                    + ax[16] * eta[i] * xi[i]
                    + ax[17] * eta[i] ** 2 * xi[i]
                    + ax[18] * eta[i] * xi[i] ** 2
                    + 2 * ax[19] * eta[i] * xi[i] * z[i],
                    ay[3]
                    + ay[5] * z[i]
                    + ay[6] * eta[i]
                    + 2 * ay[9] * z[i]
                    + ay[12] * xi[i] ** 2
                    + 2 * ay[13] * xi[i] * z[i]
                    + 2 * ay[14] * z[i] * eta[i]
                    + ay[15] * eta[i] ** 2
                    + ay[16] * eta[i] * xi[i]
                    + ay[17] * eta[i] ** 2 * xi[i]
                    + ay[18] * eta[i] * xi[i] ** 2
                    + 2 * ay[19] * eta[i] * xi[i] * z[i],
                    az[3]
                    + az[5] * z[i]
                    + az[6] * eta[i]
                    + 2 * az[9] * z[i]
                    + az[12] * xi[i] ** 2
                    + 2 * az[13] * xi[i] * z[i]
                    + 2 * az[14] * z[i] * eta[i]
                    + az[15] * eta[i] ** 2
                    + az[16] * eta[i] * xi[i]
                    + az[17] * eta[i] ** 2 * xi[i]
                    + az[18] * eta[i] * xi[i] ** 2
                    + 2 * az[19] * eta[i] * xi[i] * z[i],
                ],
            ]
        )
    return J


def _B2D4(xi: np.ndarray, eta: np.ndarray, x: np.ndarray, nInt: int):
    """Get the B operator for a Quad4 element."""

    Bi = np.zeros([nInt, 3, 8])
    # [a] matrix that connects strain and displacement derivatives
    a = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 1, 0]])
    J = _J2D4(xi, eta, x, nInt)
    for i in range(0, nInt):  # for all Gauss points (N in total)
        # make inverse of Jacobi
        invJ = np.linalg.inv(J[i])
        # [b] connects displacement derivatives (Q4)
        bi = np.array([[invJ, np.zeros([2, 2])], [np.zeros([2, 2]), invJ]])
        # make [b] what it should actually look like
        b = bi.transpose(0, 2, 1, 3).reshape(4, 4)
        # [h] as a temporary matrix
        h = np.array(
            [
                [
                    -1 / 4 * (1 - xi[i]),
                    0,
                    1 / 4 * (1 - xi[i]),
                    0,
                    1 / 4 * (1 + xi[i]),
                    0,
                    -1 / 4 * (1 + xi[i]),
                ],
                [
                    -1 / 4 * (1 - eta[i]),
                    0,
                    -1 / 4 * (1 + eta[i]),
                    0,
                    1 / 4 * (1 + eta[i]),
                    0,
                    1 / 4 * (1 - eta[i]),
                ],
            ]
        )
        # assemble [c] differentiated shapefunctions (Q4)
        c = np.vstack([np.hstack([h, np.zeros([2, 1])]), np.hstack([np.zeros([2, 1]), h])])
        # [B] for all different s and t
        Bi[i] = a @ b @ c
    return Bi


def _B2D8(xi: np.ndarray, eta: np.ndarray, x: np.ndarray, nInt: int):
    """Get the B operator for a Quad8 element."""

    Bi = np.zeros([nInt, 3, 16])
    # [a] matrix that connects strain and displacement derivatives
    a = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 1, 0]])
    J = _J2D8(xi, eta, x, nInt)
    for i in range(0, nInt):  # for all Gauss points (N in total)
        # make inverse of Jacobi
        invJ = np.linalg.inv(J[i])
        # [b] connects displacement derivatives (Q8)
        bi = np.array([[invJ, np.zeros([2, 2])], [np.zeros([2, 2]), invJ]])
        # make [b] what it should actually look like
        b = bi.transpose(0, 2, 1, 3).reshape(4, 4)
        # [h] as a temporary matrix
        h = np.array(
            [
                [
                    -1 / 4 * (-1 + xi[i]) * (2 * eta[i] + xi[i]),
                    0,
                    1 / 4 * (-1 + xi[i]) * (xi[i] - 2 * eta[i]),
                    0,
                    1 / 4 * (1 + xi[i]) * (2 * eta[i] + xi[i]),
                    0,
                    -1 / 4 * (1 + xi[i]) * (xi[i] - 2 * eta[i]),
                    0,
                    eta[i] * (-1 + xi[i]),
                    0,
                    -1 / 2 * (1 + xi[i]) * (-1 + xi[i]),
                    0,
                    -eta[i] * (1 + xi[i]),
                    0,
                    1 / 2 * (1 + xi[i]) * (-1 + xi[i]),
                ],
                [
                    -1 / 4 * (-1 + eta[i]) * (eta[i] + 2 * xi[i]),
                    0,
                    1 / 4 * (1 + eta[i]) * (2 * xi[i] - eta[i]),
                    0,
                    1 / 4 * (1 + eta[i]) * (eta[i] + 2 * xi[i]),
                    0,
                    -1 / 4 * (-1 + eta[i]) * (2 * xi[i] - eta[i]),
                    0,
                    1 / 2 * (1 + eta[i]) * (-1 + eta[i]),
                    0,
                    -xi[i] * (1 + eta[i]),
                    0,
                    -1 / 2 * (1 + eta[i]) * (-1 + eta[i]),
                    0,
                    xi[i] * (-1 + eta[i]),
                ],
            ]
        )
        # assemble [c] differentiated shapefunctions (Q8)
        c = np.vstack([np.hstack([h, np.zeros([2, 1])]), np.hstack([np.zeros([2, 1]), h])])
        # [B] for all different s and t
        Bi[i] = a @ b @ c
    return Bi


def _B3D8(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, x: np.ndarray, nInt: int):
    """Get the B operator for a Hexa8 element."""

    Bi = np.zeros([nInt, 6, 24])
    # [a] matrix that connects strain and displacement derivatives
    a = np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 1, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 1, 0],
        ]
    )
    J = _J3D8(xi, eta, z, x, nInt)
    for i in range(0, nInt):  # for all Gauss points (N in total)
        # make inverse of Jacobi
        invJ = np.linalg.inv(J[i])
        # [b] connects displacement derivatives (Hex8)
        bi = np.array(
            [
                [invJ, np.zeros([3, 3]), np.zeros([3, 3])],
                [np.zeros([3, 3]), invJ, np.zeros([3, 3])],
                [np.zeros([3, 3]), np.zeros([3, 3]), invJ],
            ]
        )
        # make [b] what it should actually look like
        b = bi.transpose(0, 2, 1, 3).reshape(9, 9)
        # [h] as a temporary matrix
        h = (
            1
            / 8
            * np.array(
                [
                    [
                        -(1 - xi[i]) * (1 - z[i]),
                        0,
                        0,
                        -(1 - xi[i]) * (1 + z[i]),
                        0,
                        0,
                        (1 - xi[i]) * (1 + z[i]),
                        0,
                        0,
                        (1 - xi[i]) * (1 - z[i]),
                        0,
                        0,
                        -(1 + xi[i]) * (1 - z[i]),
                        0,
                        0,
                        -(1 + xi[i]) * (1 + z[i]),
                        0,
                        0,
                        (1 + xi[i]) * (1 + z[i]),
                        0,
                        0,
                        (1 + xi[i]) * (1 - z[i]),
                    ],
                    [
                        -(1 - eta[i]) * (1 - z[i]),
                        0,
                        0,
                        -(1 - eta[i]) * (1 + z[i]),
                        0,
                        0,
                        -(1 + eta[i]) * (1 + z[i]),
                        0,
                        0,
                        -(1 + eta[i]) * (1 - z[i]),
                        0,
                        0,
                        (1 - eta[i]) * (1 - z[i]),
                        0,
                        0,
                        (1 - eta[i]) * (1 + z[i]),
                        0,
                        0,
                        (1 + eta[i]) * (1 + z[i]),
                        0,
                        0,
                        (1 + eta[i]) * (1 - z[i]),
                    ],
                    [
                        -(1 - eta[i]) * (1 - xi[i]),
                        0,
                        0,
                        (1 - eta[i]) * (1 - xi[i]),
                        0,
                        0,
                        (1 + eta[i]) * (1 - xi[i]),
                        0,
                        0,
                        -(1 + eta[i]) * (1 - xi[i]),
                        0,
                        0,
                        -(1 - eta[i]) * (1 + xi[i]),
                        0,
                        0,
                        (1 - eta[i]) * (1 + xi[i]),
                        0,
                        0,
                        (1 + eta[i]) * (1 + xi[i]),
                        0,
                        0,
                        -(1 + eta[i]) * (1 + xi[i]),
                    ],
                ]
            )
        )
        # assemble [c] differentiated shapefunctions (Hexa8)
        c = np.vstack(
            [
                np.hstack([h, np.zeros([3, 2])]),
                np.hstack([np.zeros([3, 1]), h, np.zeros([3, 1])]),
                np.hstack([np.zeros([3, 2]), h]),
            ]
        )
        # [B] for all different s and t
        Bi[i] = a @ b @ c
    return Bi


def _B3D20(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, x: np.ndarray, nInt: int):
    """Get the B operator for a Hexa20 element."""

    Bi = np.zeros([nInt, 6, 60])
    # [a] matrix that connects strain and displacement derivatives
    a = np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 1, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 1, 0],
        ]
    )
    J = _J3D20(xi, eta, z, x, nInt)
    for i in range(0, nInt):  # for all Gauss points (N in total)
        # make inverse of Jacobi
        invJ = np.linalg.inv(J[i])
        # make [b] - connects displacement derivatives (Hexa8)
        b = np.bmat(
            [
                [invJ, np.zeros([3, 3]), np.zeros([3, 3])],
                [np.zeros([3, 3]), invJ, np.zeros([3, 3])],
                [np.zeros([3, 3]), np.zeros([3, 3]), invJ],
            ]
        )
        # [h] as a temporary matrix
        h = np.array(
            [
                [
                    ((xi[i] - 1) * (z[i] - 1) * (2 * eta[i] + xi[i] + z[i] + 1)) / 8,
                    0,
                    0,
                    -((xi[i] - 1) * (z[i] + 1) * (2 * eta[i] + xi[i] - z[i] + 1)) / 8,
                    0,
                    0,
                    -((xi[i] - 1) * (z[i] + 1) * (2 * eta[i] - xi[i] + z[i] - 1)) / 8,
                    0,
                    0,
                    -((xi[i] - 1) * (z[i] - 1) * (xi[i] - 2 * eta[i] + z[i] + 1)) / 8,
                    0,
                    0,
                    -((xi[i] + 1) * (z[i] - 1) * (2 * eta[i] - xi[i] + z[i] + 1)) / 8,
                    0,
                    0,
                    ((xi[i] + 1) * (z[i] + 1) * (2 * eta[i] - xi[i] - z[i] + 1)) / 8,
                    0,
                    0,
                    ((xi[i] + 1) * (z[i] + 1) * (2 * eta[i] + xi[i] + z[i] - 1)) / 8,
                    0,
                    0,
                    -((xi[i] + 1) * (z[i] - 1) * (2 * eta[i] + xi[i] - z[i] - 1)) / 8,
                    0,
                    0,
                    -((z[i] ** 2 - 1) * (xi[i] - 1)) / 4,
                    0,
                    0,
                    (eta[i] * (xi[i] - 1) * (z[i] + 1)) / 2,
                    0,
                    0,
                    ((z[i] ** 2 - 1) * (xi[i] - 1)) / 4,
                    0,
                    0,
                    -(eta[i] * (xi[i] - 1) * (z[i] - 1)) / 2,
                    0,
                    0,
                    ((z[i] ** 2 - 1) * (xi[i] + 1)) / 4,
                    0,
                    0,
                    -(eta[i] * (xi[i] + 1) * (z[i] + 1)) / 2,
                    0,
                    0,
                    -((z[i] ** 2 - 1) * (xi[i] + 1)) / 4,
                    0,
                    0,
                    (eta[i] * (xi[i] + 1) * (z[i] - 1)) / 2,
                    0,
                    0,
                    -((xi[i] ** 2 - 1) * (z[i] - 1)) / 4,
                    0,
                    0,
                    ((xi[i] ** 2 - 1) * (z[i] + 1)) / 4,
                    0,
                    0,
                    -((xi[i] ** 2 - 1) * (z[i] + 1)) / 4,
                    0,
                    0,
                    ((xi[i] ** 2 - 1) * (z[i] - 1)) / 4,
                ],
                [
                    ((eta[i] - 1) * (z[i] - 1) * (eta[i] + 2 * xi[i] + z[i] + 1)) / 8,
                    0,
                    0,
                    -((eta[i] - 1) * (z[i] + 1) * (eta[i] + 2 * xi[i] - z[i] + 1)) / 8,
                    0,
                    0,
                    -((eta[i] + 1) * (z[i] + 1) * (eta[i] - 2 * xi[i] + z[i] - 1)) / 8,
                    0,
                    0,
                    -((eta[i] + 1) * (z[i] - 1) * (2 * xi[i] - eta[i] + z[i] + 1)) / 8,
                    0,
                    0,
                    -((eta[i] - 1) * (z[i] - 1) * (eta[i] - 2 * xi[i] + z[i] + 1)) / 8,
                    0,
                    0,
                    ((eta[i] - 1) * (z[i] + 1) * (eta[i] - 2 * xi[i] - z[i] + 1)) / 8,
                    0,
                    0,
                    ((eta[i] + 1) * (z[i] + 1) * (eta[i] + 2 * xi[i] + z[i] - 1)) / 8,
                    0,
                    0,
                    -((eta[i] + 1) * (z[i] - 1) * (eta[i] + 2 * xi[i] - z[i] - 1)) / 8,
                    0,
                    0,
                    -(eta[i] / 4 - 1 / 4) * (z[i] ** 2 - 1),
                    0,
                    0,
                    (eta[i] ** 2 / 4 - 1 / 4) * (z[i] + 1),
                    0,
                    0,
                    (eta[i] / 4 + 1 / 4) * (z[i] ** 2 - 1),
                    0,
                    0,
                    -(eta[i] ** 2 / 4 - 1 / 4) * (z[i] - 1),
                    0,
                    0,
                    (eta[i] / 4 - 1 / 4) * (z[i] ** 2 - 1),
                    0,
                    0,
                    -(eta[i] ** 2 / 4 - 1 / 4) * (z[i] + 1),
                    0,
                    0,
                    -(eta[i] / 4 + 1 / 4) * (z[i] ** 2 - 1),
                    0,
                    0,
                    (eta[i] ** 2 / 4 - 1 / 4) * (z[i] - 1),
                    0,
                    0,
                    -2 * xi[i] * (eta[i] / 4 - 1 / 4) * (z[i] - 1),
                    0,
                    0,
                    2 * xi[i] * (eta[i] / 4 - 1 / 4) * (z[i] + 1),
                    0,
                    0,
                    -2 * xi[i] * (eta[i] / 4 + 1 / 4) * (z[i] + 1),
                    0,
                    0,
                    2 * xi[i] * (eta[i] / 4 + 1 / 4) * (z[i] - 1),
                ],
                [
                    ((eta[i] - 1) * (xi[i] - 1) * (eta[i] + xi[i] + 2 * z[i] + 1)) / 8,
                    0,
                    0,
                    -((eta[i] - 1) * (xi[i] - 1) * (eta[i] + xi[i] - 2 * z[i] + 1)) / 8,
                    0,
                    0,
                    -((eta[i] + 1) * (xi[i] - 1) * (eta[i] - xi[i] + 2 * z[i] - 1)) / 8,
                    0,
                    0,
                    -((eta[i] + 1) * (xi[i] - 1) * (xi[i] - eta[i] + 2 * z[i] + 1)) / 8,
                    0,
                    0,
                    -((eta[i] - 1) * (xi[i] + 1) * (eta[i] - xi[i] + 2 * z[i] + 1)) / 8,
                    0,
                    0,
                    ((eta[i] - 1) * (xi[i] + 1) * (eta[i] - xi[i] - 2 * z[i] + 1)) / 8,
                    0,
                    0,
                    ((eta[i] + 1) * (xi[i] + 1) * (eta[i] + xi[i] + 2 * z[i] - 1)) / 8,
                    0,
                    0,
                    -((eta[i] + 1) * (xi[i] + 1) * (eta[i] + xi[i] - 2 * z[i] - 1)) / 8,
                    0,
                    0,
                    -2 * z[i] * (eta[i] / 4 - 1 / 4) * (xi[i] - 1),
                    0,
                    0,
                    (eta[i] ** 2 / 4 - 1 / 4) * (xi[i] - 1),
                    0,
                    0,
                    2 * z[i] * (eta[i] / 4 + 1 / 4) * (xi[i] - 1),
                    0,
                    0,
                    -(eta[i] ** 2 / 4 - 1 / 4) * (xi[i] - 1),
                    0,
                    0,
                    2 * z[i] * (eta[i] / 4 - 1 / 4) * (xi[i] + 1),
                    0,
                    0,
                    -(eta[i] ** 2 / 4 - 1 / 4) * (xi[i] + 1),
                    0,
                    0,
                    -2 * z[i] * (eta[i] / 4 + 1 / 4) * (xi[i] + 1),
                    0,
                    0,
                    (eta[i] ** 2 / 4 - 1 / 4) * (xi[i] + 1),
                    0,
                    0,
                    -(eta[i] / 4 - 1 / 4) * (xi[i] ** 2 - 1),
                    0,
                    0,
                    (eta[i] / 4 - 1 / 4) * (xi[i] ** 2 - 1),
                    0,
                    0,
                    -(eta[i] / 4 + 1 / 4) * (xi[i] ** 2 - 1),
                    0,
                    0,
                    (eta[i] / 4 + 1 / 4) * (xi[i] ** 2 - 1),
                ],
            ]
        )
        # assemble [c] differentiated shapefunctions (Hexa20)
        c = np.vstack(
            [
                np.hstack([h, np.zeros([3, 2])]),
                np.hstack([np.zeros([3, 1]), h, np.zeros([3, 1])]),
                np.hstack([np.zeros([3, 2]), h]),
            ]
        )
        # [B] for all different s and t
        Bi[i] = a @ b @ c
    return Bi
