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
import numpy.linalg as lin


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
        Dimension the element has.

    Returns
    -------
    np.ndarray
        The requested Jacobian matrix."""

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


def computeDeformationGradient(U: np.ndarray, nablaN: np.ndarray, nInt: int, nNodes: int, dim: int):
    """Get the deformation gradient for a nonlinear element.

    Parameters
    ----------
    U
        The current displacement vector.
    nablaN
        The derivative of the shape functions w.r.t the actual coordinates.
    nInt
        Number of quadrature points the element has.
    nNodes
        Number of nodes the element has.
    dim
        Dimension of the domain.

    Returns
    -------
    np.ndarray
        The deformation gradient."""

    F = np.zeros([nInt, dim, dim])
    if dim == 2:
        for i in range(nInt):  # for all Gauss points (N in total)
            H = np.zeros([2, 2])
            for j in range(nNodes):  # make deformation gradient
                UI = U[2 * j : 2 * (j + 1)]
                H += np.outer(UI, nablaN[i, :, j])
            F[i] = np.eye(2) + H
    elif dim == 3:
        for i in range(nInt):  # for all Gauss points (N in total)
            H = np.zeros([3, 3])
            for j in range(nNodes):  # for all points/shape functions
                UI = U[3 * j : 3 * (j + 1)]
                H += np.outer(UI, nablaN[i, :, j])
            F[i] = np.eye(3) + H
    return F


def computeBOperator(F: np.ndarray, nablaN: np.ndarray, nInt: int, nNodes: int, dim: int):
    """Get the B operator for the element calculation.

    Parameters
    ----------
    F
        The deformation gradient.
    nablaN
        The derivative of the shape functions w.r.t the actual coordinates.
    nInt
        Number of integration points.
    nNodes
        Number of nodes the element has.
    dim
        Dimension the element has.

    Returns
    -------
    np.ndarray
        The requested B operator for a nonlinear element.
    np.ndarray
        The deformation gradient."""

    if dim == 2:
        return _B02D(F, nablaN, nInt, nNodes)
    elif dim == 3:
        return _B03D(F, nablaN, nInt, nNodes)


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

    Returns
    -------
    np.ndarray
        The shape functions at the given coordinates."""

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


def computeNablaN(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, J: np.ndarray, nInt: int, nNodes: int, dim: int):
    """Get the nabla N(i) = dN/dX operator.

    Parameters
    ----------
    xi
        Local coordinates xi for the integration points.
    eta
        Local coordinates eta for the integration points.
    z
        Local coordinates zeta for the integration points.
    J
        The jacobian between the local and global undeformed coordinates.
    nInt
        Number of integration points.
    nNodes
        Number of nodes the element has.
    dim
        Dimension the element has.

    Returns
    -------
    np.ndarray
        The derivative of the shape functions w.r.t the actual coordinates."""

    nablaN = np.zeros([nInt, dim, nNodes])
    for i in range(nInt):  # for all Gauss points (N in total)
        if dim == 2:
            dN = _Ndiff2D(xi[i], eta[i], nNodes)
        else:
            dN = _Ndiff3D(xi[i], eta[i], z[i], nNodes)
        invJ = lin.inv(J[i])
        for j in range(nNodes):  # for all points/shape functions
            nablaN[i, :, j] = dN[:, j] @ invJ.T
    return nablaN


def makeH2D(Hmat, dim):
    """Puts the element material stiffness from 4D into (nDof, nDof).

    Parameters
    ----------
    Hmat
        The element material stiffness matrix.
    dim
        Dimension of the domain.

    Returns
    -------
    np.ndarray
        The element stiffness matrix in 2D form."""

    if dim == 2:
        Hmat = Hmat[:, :dim, :, :dim]
    m, n, _, _ = Hmat.shape
    HmatDim = np.zeros([m * n, m * n])
    for i in range(m):
        for j in range(n):
            HmatDim[dim * i + j] = Hmat[i, j].flatten()
    return HmatDim


def _Ndiff2D(xi: np.ndarray, eta: np.ndarray, nNodes: int):
    """Calculate the differentiated 2D shape functions.

    Parameters
    ----------
    xi
        Local coordinates xi for the integration points.
    eta
        Local coordinates eta for the integration points.
    nNodes
        Number of nodes the element has.

    Returns
    -------
    np.ndarray
        The derivative of the 2D shape functions w.r.t. the local coordinates."""

    if nNodes == 4:  # Quad4
        return np.array(
            [
                [
                    -1 / 4 * (1 - xi),
                    1 / 4 * (1 - xi),
                    1 / 4 * (1 + xi),
                    -1 / 4 * (1 + xi),
                ],
                [
                    -1 / 4 * (1 - eta),
                    -1 / 4 * (1 + eta),
                    1 / 4 * (1 + eta),
                    1 / 4 * (1 - eta),
                ],
            ]
        )
    else:  # Quad8
        return np.array(
            [
                [
                    -1 / 4 * (-1 + xi) * (2 * eta + xi),
                    1 / 4 * (-1 + xi) * (xi - 2 * eta),
                    1 / 4 * (1 + xi) * (2 * eta + xi),
                    -1 / 4 * (1 + xi) * (xi - 2 * eta),
                    eta * (-1 + xi),
                    -1 / 2 * (1 + xi) * (-1 + xi),
                    -eta * (1 + xi),
                    1 / 2 * (1 + xi) * (-1 + xi),
                ],
                [
                    -1 / 4 * (-1 + eta) * (eta + 2 * xi),
                    1 / 4 * (1 + eta) * (2 * xi - eta),
                    1 / 4 * (1 + eta) * (eta + 2 * xi),
                    -1 / 4 * (-1 + eta) * (2 * xi - eta),
                    1 / 2 * (1 + eta) * (-1 + eta),
                    -xi * (1 + eta),
                    -1 / 2 * (1 + eta) * (-1 + eta),
                    xi * (-1 + eta),
                ],
            ]
        )


def _Ndiff3D(xi: np.ndarray, eta: np.ndarray, z: np.ndarray, nNodes: int):
    """Calculate the differentiated 3D shape functions.

    Parameters
    ----------
    xi
        Local coordinates xi for the integration points.
    eta
        Local coordinates eta for the integration points.
    z
        Local coordinates zeta for the integration points.
    nNodes
        Number of nodes the element has.

    Returns
    -------
    np.ndarray
        The derivative of the 3D shape functions w.r.t. the local coordinates."""

    if nNodes == 8:  # Hexa8
        return (
            1
            / 8
            * np.array(
                [
                    [
                        -(1 - xi) * (1 - z),
                        -(1 - xi) * (1 + z),
                        (1 - xi) * (1 + z),
                        (1 - xi) * (1 - z),
                        -(1 + xi) * (1 - z),
                        -(1 + xi) * (1 + z),
                        (1 + xi) * (1 + z),
                        (1 + xi) * (1 - z),
                    ],
                    [
                        -(1 - eta) * (1 - z),
                        -(1 - eta) * (1 + z),
                        -(1 + eta) * (1 + z),
                        -(1 + eta) * (1 - z),
                        (1 - eta) * (1 - z),
                        (1 - eta) * (1 + z),
                        (1 + eta) * (1 + z),
                        (1 + eta) * (1 - z),
                    ],
                    [
                        -(1 - eta) * (1 - xi),
                        (1 - eta) * (1 - xi),
                        (1 + eta) * (1 - xi),
                        -(1 + eta) * (1 - xi),
                        -(1 - eta) * (1 + xi),
                        (1 - eta) * (1 + xi),
                        (1 + eta) * (1 + xi),
                        -(1 + eta) * (1 + xi),
                    ],
                ]
            )
        )
    else:  # Hexa20
        return np.array(
            [
                [
                    ((xi - 1) * (z - 1) * (2 * eta + xi + z + 1)) / 8,
                    -((xi - 1) * (z + 1) * (2 * eta + xi - z + 1)) / 8,
                    -((xi - 1) * (z + 1) * (2 * eta - xi + z - 1)) / 8,
                    -((xi - 1) * (z - 1) * (xi - 2 * eta + z + 1)) / 8,
                    -((xi + 1) * (z - 1) * (2 * eta - xi + z + 1)) / 8,
                    ((xi + 1) * (z + 1) * (2 * eta - xi - z + 1)) / 8,
                    ((xi + 1) * (z + 1) * (2 * eta + xi + z - 1)) / 8,
                    -((xi + 1) * (z - 1) * (2 * eta + xi - z - 1)) / 8,
                    -((z**2 - 1) * (xi - 1)) / 4,
                    (eta * (xi - 1) * (z + 1)) / 2,
                    ((z**2 - 1) * (xi - 1)) / 4,
                    -(eta * (xi - 1) * (z - 1)) / 2,
                    ((z**2 - 1) * (xi + 1)) / 4,
                    -(eta * (xi + 1) * (z + 1)) / 2,
                    -((z**2 - 1) * (xi + 1)) / 4,
                    (eta * (xi + 1) * (z - 1)) / 2,
                    -((xi**2 - 1) * (z - 1)) / 4,
                    ((xi**2 - 1) * (z + 1)) / 4,
                    -((xi**2 - 1) * (z + 1)) / 4,
                    ((xi**2 - 1) * (z - 1)) / 4,
                ],
                [
                    ((eta - 1) * (z - 1) * (eta + 2 * xi + z + 1)) / 8,
                    -((eta - 1) * (z + 1) * (eta + 2 * xi - z + 1)) / 8,
                    -((eta + 1) * (z + 1) * (eta - 2 * xi + z - 1)) / 8,
                    -((eta + 1) * (z - 1) * (2 * xi - eta + z + 1)) / 8,
                    -((eta - 1) * (z - 1) * (eta - 2 * xi + z + 1)) / 8,
                    ((eta - 1) * (z + 1) * (eta - 2 * xi - z + 1)) / 8,
                    ((eta + 1) * (z + 1) * (eta + 2 * xi + z - 1)) / 8,
                    -((eta + 1) * (z - 1) * (eta + 2 * xi - z - 1)) / 8,
                    -(eta / 4 - 1 / 4) * (z**2 - 1),
                    (eta**2 / 4 - 1 / 4) * (z + 1),
                    (eta / 4 + 1 / 4) * (z**2 - 1),
                    -(eta**2 / 4 - 1 / 4) * (z - 1),
                    (eta / 4 - 1 / 4) * (z**2 - 1),
                    -(eta**2 / 4 - 1 / 4) * (z + 1),
                    -(eta / 4 + 1 / 4) * (z**2 - 1),
                    (eta**2 / 4 - 1 / 4) * (z - 1),
                    -2 * xi * (eta / 4 - 1 / 4) * (z - 1),
                    2 * xi * (eta / 4 - 1 / 4) * (z + 1),
                    -2 * xi * (eta / 4 + 1 / 4) * (z + 1),
                    2 * xi * (eta / 4 + 1 / 4) * (z - 1),
                ],
                [
                    ((eta - 1) * (xi - 1) * (eta + xi + 2 * z + 1)) / 8,
                    -((eta - 1) * (xi - 1) * (eta + xi - 2 * z + 1)) / 8,
                    -((eta + 1) * (xi - 1) * (eta - xi + 2 * z - 1)) / 8,
                    -((eta + 1) * (xi - 1) * (xi - eta + 2 * z + 1)) / 8,
                    -((eta - 1) * (xi + 1) * (eta - xi + 2 * z + 1)) / 8,
                    ((eta - 1) * (xi + 1) * (eta - xi - 2 * z + 1)) / 8,
                    ((eta + 1) * (xi + 1) * (eta + xi + 2 * z - 1)) / 8,
                    -((eta + 1) * (xi + 1) * (eta + xi - 2 * z - 1)) / 8,
                    -2 * z * (eta / 4 - 1 / 4) * (xi - 1),
                    (eta**2 / 4 - 1 / 4) * (xi - 1),
                    2 * z * (eta / 4 + 1 / 4) * (xi - 1),
                    -(eta**2 / 4 - 1 / 4) * (xi - 1),
                    2 * z * (eta / 4 - 1 / 4) * (xi + 1),
                    -(eta**2 / 4 - 1 / 4) * (xi + 1),
                    -2 * z * (eta / 4 + 1 / 4) * (xi + 1),
                    (eta**2 / 4 - 1 / 4) * (xi + 1),
                    -(eta / 4 - 1 / 4) * (xi**2 - 1),
                    (eta / 4 - 1 / 4) * (xi**2 - 1),
                    -(eta / 4 + 1 / 4) * (xi**2 - 1),
                    (eta / 4 + 1 / 4) * (xi**2 - 1),
                ],
            ]
        )


def _J2D4(xi: np.ndarray, eta: np.ndarray, x: np.ndarray, nInt: int):
    """Get the Jacobi matrix for a Quad4 element.

    Parameters
    ----------
    xi
        Local coordinates xi for the integration points.
    eta
        Local coordinates eta for the integration points.
    x
        Coordinates of the element points.
    nInt
        Number of quadrature points the element has.

    Returns
    -------
    np.ndarray
        The requested Jacobian matrix."""

    J = np.zeros([nInt, 2, 2])
    # calc all parameters for the X and Y functions (Q4)
    A = np.array([[1, -1, -1, 1], [1, 1, -1, -1], [1, 1, 1, 1], [1, -1, 1, -1]])
    invA = lin.solve(A, np.eye(4))
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
    """Get the Jacobi matrix for a Quad8 element.

    Parameters
    ----------
    xi
        Local coordinates xi for the integration points.
    eta
        Local coordinates eta for the integration points.
    x
        Coordinates of the element points.
    nInt
        Number of quadrature points the element has.

    Returns
    -------
    np.ndarray
        The requested Jacobian matrix."""

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
    invA = lin.solve(A, np.eye(8))
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
    """Get the Jacobi matrix for a Hexa8 element.

    Parameters
    ----------
    xi
        Local coordinates xi for the integration points.
    eta
        Local coordinates eta for the integration points.
    z
        Local coordinates zeta for the integration points.
    x
        Coordinates of the element points.
    nInt
        Number of quadrature points the element has.

    Returns
    -------
    np.ndarray
        The requested Jacobian matrix."""

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
    invA = lin.solve(A, np.eye(8))
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
    """Get the Jacobi matrix for a Hexa20 element.

    Parameters
    ----------
    xi
        Local coordinates xi for the integration points.
    eta
        Local coordinates eta for the integration points.
    z
        Local coordinates zeta for the integration points.
    x
        Coordinates of the element points.
    nInt
        Number of quadrature points the element has.

    Returns
    -------
    np.ndarray
        The requested Jacobian matrix."""

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
    invA = lin.solve(A, np.eye(20))
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


# B0 operator
def _B02D(F: np.ndarray, nablaN: np.ndarray, nInt: int, nNodes: int):
    """Get the B operator and the deformation gradient for a nonlinear Quad element.

    Parameters
    ----------
    F
        The deformation gradient.
    nablaN
        The derivative of the shape functions w.r.t the actual coordinates.
    nInt
        Number of quadrature points the element has.
    nNodes
        Number of nodes the element has.

    Returns
    -------
    np.ndarray
        The requested B operator for a nonlinear element.
    np.ndarray
        The deformation gradient."""

    Bi = np.zeros([nInt, 3, nNodes * 2])
    for i in range(nInt):  # for all Gauss points (N in total)
        # [B] for all different xi and eta
        dxdX = F[i]
        for j in range(nNodes):  # construct B operator
            dNdX, dNdY = nablaN[i, :, j]
            Bi[i, :, 2 * j : 2 * j + 2] = np.array(
                [
                    [dNdX * dxdX[0, 0], dNdX * dxdX[1, 0]],
                    [dNdY * dxdX[0, 1], dNdY * dxdX[1, 1]],
                    [dNdX * dxdX[0, 1] + dNdY * dxdX[0, 0], dNdX * dxdX[1, 1] + dNdY * dxdX[1, 0]],
                ]
            )
    return Bi


def _B03D(F: np.ndarray, nablaN: np.ndarray, nInt: int, nNodes: int):
    """Get the B operator and the deformation gradient for a nonlinear Hexa element.

    Parameters
    ----------
    F
        The deformation gradient.
    nablaN
        The derivative of the shape functions w.r.t the actual coordinates.
    nInt
        Number of quadrature points the element has.
    nNodes
        Number of nodes the element has.

    Returns
    -------
    np.ndarray
        The requested B operator for a nonlinear element.
    np.ndarray
        The deformation gradient."""

    Bi = np.zeros([nInt, 6, nNodes * 3])
    for i in range(nInt):  # for all Gauss points (N in total)
        # [B] for all different xi, eta and zeta
        dxdX = F[i]
        for j in range(nNodes):  # construct B operator
            dNdX, dNdY, dNdZ = nablaN[i, :, j]
            Bi[i, :, 3 * j : 3 * j + 3] = np.array(
                [
                    [dNdX * dxdX[0, 0], dNdX * dxdX[1, 0], dNdX * dxdX[2, 0]],
                    [dNdY * dxdX[0, 1], dNdY * dxdX[1, 1], dNdY * dxdX[2, 1]],
                    [dNdZ * dxdX[0, 2], dNdZ * dxdX[1, 2], dNdZ * dxdX[2, 2]],
                    [
                        dNdX * dxdX[0, 1] + dNdY * dxdX[0, 0],
                        dNdX * dxdX[1, 1] + dNdY * dxdX[1, 0],
                        dNdX * dxdX[2, 1] + dNdY * dxdX[2, 0],
                    ],
                    [
                        dNdY * dxdX[0, 2] + dNdZ * dxdX[0, 1],
                        dNdY * dxdX[1, 2] + dNdZ * dxdX[1, 1],
                        dNdY * dxdX[2, 2] + dNdZ * dxdX[2, 1],
                    ],
                    [
                        dNdX * dxdX[0, 2] + dNdZ * dxdX[0, 0],
                        dNdX * dxdX[1, 2] + dNdZ * dxdX[1, 0],
                        dNdX * dxdX[2, 2] + dNdZ * dxdX[2, 0],
                    ],
                ]
            )
    return Bi
