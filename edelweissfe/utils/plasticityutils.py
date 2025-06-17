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

from edelweissfe.utils.voigtnotation import doVoigtStress

tol = 1e-16


def deviatoric(stress):
    """Computes the deviatoric stress.

    Parameters
    ----------
    stress
        The given stress as a matrix.

    Returns
    -------
    np.ndarray
        The deviatoric stress."""

    return stress - np.trace(stress) / 3 * np.eye(3)


def deviatoricNorm(devStress):
    """Computes the norm of the deviatoric stress.

    Parameters
    ----------
    devStress
        The given deviatoric stress as a matrix.

    Returns
    -------
    np.ndarray
        The norm of the deviatoric stress."""

    devStress = doVoigtStress(3, devStress)  # put in Voigt
    return np.sqrt(np.sum(np.square(devStress[0:3])) + 2 * np.sum(np.square(devStress[3:6])))


def tensorLogarithmEig(A):
    """Computes the logarithm of a tensor using eigenvalues and eigenvectors.

    Parameters
    ----------
    A
        The tensor, the logarithm should be computed of.

    Returns
    -------
    np.ndarray
        The tensor logarithm of A."""

    (x, V) = lin.eig(A)
    A_ = np.diag(np.log(x))
    return V @ A_ @ lin.inv(V)


def tensorExp(A):
    """Computes the tensor exponential using an iterative method.

    Parameters
    ----------
    A
        The tensor, the exponential should be computed of.

    Returns
    -------
    np.ndarray
        The tensor exponential of A.
    int
        The number of additive terms of the tensor exponential."""

    n = 0
    expX = np.eye(3)
    Xn = A
    nF = np.prod(np.arange(1, n + 1))
    while lin.norm(Xn) / nF > tol:
        n += 1
        nF = np.prod(np.arange(1, n + 1))
        Xn = lin.matrix_power(A, n)
        expX += Xn / nF
    return (expX, n)


def dTensorExp_dA(A, n):
    """Computes the derivative of the tensor exponential using an iterative method.

    Parameters
    ----------
    A
        The tensor at which the derivative of the tensor exponential should be computed.
    n
        The number of additive terms of the tensor exponential.

    Returns
    -------
    np.ndarray
        The 4D derivative of the tensor exponential at A."""

    dExpX = np.zeros([3, 3, 3, 3])
    for i in range(1, n + 1):
        nF = np.prod(np.arange(1, i + 1))
        for m in range(1, i + 1):
            Xm1 = lin.matrix_power(A, m - 1)
            Xnm = lin.matrix_power(A, i - m)
            dExpX += 1 / nF * np.einsum("ik,lj->jilk", Xm1, Xnm)
    return dExpX


def computeEigenprojections(A):
    """Computes an eigenprojection of A.

    Parameters
    ----------
    A
        The tensor, the eigenprojection will be computed of.
    eig
        The eigenvalue affiliated with the eigenprojection.

    Returns
    -------
    np.ndarray
        All eigenprojections of A.
    np.ndarray
        All eigenvalues of A."""

    I1 = np.trace(A)
    I2 = 1 / 2 * (I1**2 - np.trace(A @ A))
    I3 = lin.det(A)

    def eigenprojection(eig):
        return eig / (2 * eig**3 - I1 * eig**2 + I3) * (A @ A - (I1 - eig) * A + I3 / eig * np.eye(3))

    # eigenvalue computation
    R = (-2 * I1**3 + 9 * I1 * I2 - 27 * I3) / 54
    Q = np.abs(I1**2 - 3 * I2) / 9
    if Q == 0:
        eig = np.ones([3]) * I1 / 3
    else:
        S = R / np.sqrt(Q**3)
        theta = np.arccos(S if not np.abs(S) > 1.0 else np.sign(S) * 1.0)
        eig = (
            -2
            * np.sqrt(Q)
            * np.array([np.cos(theta / 3), np.cos((theta + 2 * np.pi) / 3), np.cos((theta - 2 * np.pi) / 3)])
            + I1 / 3
        )
    counts = np.sum(eig == eig[0]) + np.sum(eig == eig[1]) + np.sum(eig == eig[2])
    if counts == 3:  # no eigenvalues are the same
        E = np.zeros([3, 3, 3])
        for i in range(3):  # construction of the eigenprojections
            E[i] = eigenprojection(eig[i])
    elif eig[0] == eig[1] and eig[0] != eig[2]:  # 3rd value different
        counts = 2
        E = eigenprojection(eig[2])
    elif eig[0] == eig[2] and eig[0] != eig[1]:  # 2nd value different
        counts = 1
        E = eigenprojection(eig[1])
    elif eig[1] == eig[2] and eig[1] != eig[0]:  # 1st value different
        counts = 0
        E = eigenprojection(eig[0])
    elif counts == 9:  # all values the same
        E = np.eye(3)
    return (E, eig, counts)


def tensorLogarithm(A):
    """Computes the logarithm of a tensor using eigenvalues and eigenprojections.

    Parameters
    ----------
    A
        The tensor, the logarithm should be computed of.

    Returns
    -------
    np.ndarray
        The tensor logarithm of A."""

    (E, x, counts) = computeEigenprojections(A)
    logXi = np.log(x)
    if counts == 3:  # no eigenvalues are the same
        logA = np.zeros([3, 3])
        for i in range(3):  # construction of the eigenprojections
            logA += logXi[i] * E[i]
        return logA
    elif counts in [0, 1, 2]:
        b = 0 if counts != 0 else 1
        return logXi[counts] * E + logXi[b] * (np.eye(3) - E)
    elif counts == 9:  # all values the same
        return logXi[0] * np.eye(3)


def dTensorLog_dA(A):
    """Computes the derivative of the tensor logarithm.

    Parameters
    ----------
    A
        The matrix at which the derivative of the tensor logarithm needs to be evaluated at.

    Returns
    -------
    np.ndarray
        The derivative of the tensor logarithm at A."""

    (E, x, counts) = computeEigenprojections(A)
    y = np.log(x)
    Ie = np.eye(3)
    IS = 1 / 2 * (np.einsum("ik,jl->ijkl", Ie, Ie) + np.einsum("il,jk->ijkl", Ie, Ie))
    dX2dX = (
        1
        / 2
        * (
            np.einsum("ik,lj->ijkl", Ie, A)
            + np.einsum("il,kj->ijkl", Ie, A)
            + np.einsum("jl,ik->ijkl", Ie, A)
            + np.einsum("kj,il->ijkl", Ie, A)
        )
    )
    if counts == 3:  # no eigenvalues are the same
        DlogA = np.zeros([3, 3, 3, 3])
        for i in range(3):
            # make sure a, b, c are always different
            b = 0 if (i == 1 or i == 2) else 1
            c = 2 if (i == 0 or i == 1) else 1
            # compute tensor products
            EaEa = np.tensordot(E[i], E[i], axes=0)
            EbEb = np.tensordot(E[b], E[b], axes=0)
            EcEc = np.tensordot(E[c], E[c], axes=0)
            # sum up differential
            DlogA += (
                y[i]
                / ((x[i] - x[b]) * (x[i] - x[c]))
                * (dX2dX - (x[b] + x[c]) * IS - ((x[i] - x[b]) + (x[i] - x[c])) * EaEa - (x[b] - x[c]) * (EbEb - EcEc))
                + 1 / x[i] * EaEa
            )
        return DlogA
    elif counts in [0, 1, 2]:  # two eigenvalues are the same
        b = 0 if counts != 0 else 1

        def computeS(ya, yc, xa, xc):
            # compute the s values
            s1 = (ya - yc) / (xa - xc) ** 2 - 1 / (xc * (xa - xc))
            s2 = 2 * xc * (ya - yc) / (xa - xc) ** 2 - (xa + xc) / (xa - xc) * 1 / xc
            s3 = 2 * (ya - yc) / (xa - xc) ** 3 - (1 / xa + 1 / xc) / (xa - xc) ** 2
            s4 = xc * s3
            s6 = xc**2 * s3
            return (s1, s2, s3, s4, s4, s6)

        def computeDlogA(s1, s2, s3, s4, s5, s6):
            # compute the derivative for two same values
            XX = np.tensordot(A, A, axes=0)
            XI = np.tensordot(A, Ie, axes=0)
            IX = np.tensordot(Ie, A, axes=0)
            II = np.tensordot(Ie, Ie, axes=0)
            return s1 * dX2dX - s2 * IS - s3 * XX + s4 * XI + s5 * IX - s6 * II

        allS = computeS(y[counts], y[b], x[counts], x[b])
        return computeDlogA(*allS)
    elif counts == 9:  # all eigenvalues are the same
        return 1 / x[0] * IS
