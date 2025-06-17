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


def doVoigtStrain(dim, e):
    """Put strain symmetric values into voigt form.

    Parameters
    ----------
    dim
        Dimension of the model.
    e
        Strain in matrix form to be put into voigt form.

    Returns
    -------
    np.ndarray
        The strain in Voigt notation for the corresponding dimension."""

    if dim == 2:
        return np.array([e[0, 0], e[1, 1], e[0, 1] * 2])
    elif dim == 3:
        return np.array([e[0, 0], e[1, 1], e[2, 2], e[0, 1] * 2, e[1, 2] * 2, e[0, 2] * 2])


def undoVoigtStrain(dim, e):
    """Put strain symmetric values into matrix form.

    Parameters
    ----------
    dim
        Dimension of the model.
    e
        Strain in voigt form to be put into matrix form.

    Returns
    -------
    np.ndarray
        The strain in matrix notation for the corresponding dimension."""

    if dim == 2:
        return np.array([[e[0], e[3] / 2], [e[3] / 2, e[1]]])
    elif dim == 3:
        return np.array([[e[0], e[3] / 2, e[5] / 2], [e[3] / 2, e[1], e[4] / 2], [e[5] / 2, e[4] / 2, e[2]]])


def doVoigtStress(dim, s):
    """Put stress symmetric values into voigt form.

    Parameters
    ----------
    dim
        Dimension of the model.
    s
        Stress in matrix form to be put into voigt form.

    Returns
    -------
    np.ndarray
        The stress in Voigt notation for the corresponding dimension."""

    if dim == 2:
        return np.array([s[0, 0], s[1, 1], s[0, 1]])
    elif dim == 3:
        return np.array([s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[1, 2], s[0, 2]])


def undoVoigtStress(dim, s):
    """Put stress symmetric values into matrix form.

    Parameters
    ----------
    dim
        Dimension of the model.
    s
        Stress in voigt form to be put into matrix form.

    Returns
    -------
    np.ndarray
        The strain in matrix notation for the corresponding dimension."""

    if dim == 2:
        return np.array([[s[0], s[2]], [s[2], s[1]]])
    elif dim == 3:
        return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


def doVoigtMatrix(dim, C):
    """Put the material tangent moduli into voigt form.

    Parameters
    ----------
    dim
        Dimension of the model.
    C
        4D-Matrix to be put into 2D-matrix.

    Returns
    -------
    np.ndarray
        The tangent modulus in Voigt notation for the corresponding dimension."""

    if dim == 2:
        return np.array(
            [
                [C[0, 0, 0, 0], C[0, 0, 1, 1], C[0, 0, 0, 1]],
                [C[1, 1, 0, 0], C[1, 1, 1, 1], C[1, 1, 0, 1]],
                [C[0, 1, 0, 0], C[0, 1, 1, 1], C[0, 1, 0, 1]],
            ]
        )
    elif dim == 3:
        return np.array(
            [
                [C[0, 0, 0, 0], C[0, 0, 1, 1], C[0, 0, 2, 2], C[0, 0, 0, 1], C[0, 0, 1, 2], C[0, 0, 0, 2]],
                [C[1, 1, 0, 0], C[1, 1, 1, 1], C[1, 1, 2, 2], C[1, 1, 0, 1], C[1, 1, 1, 2], C[1, 1, 0, 2]],
                [C[2, 2, 0, 0], C[2, 2, 1, 1], C[2, 2, 2, 2], C[2, 2, 0, 1], C[2, 2, 1, 2], C[2, 2, 0, 2]],
                [C[0, 1, 0, 0], C[0, 1, 1, 1], C[0, 1, 2, 2], C[0, 1, 0, 1], C[0, 1, 1, 2], C[0, 1, 0, 2]],
                [C[1, 2, 0, 0], C[1, 2, 1, 1], C[1, 2, 2, 2], C[1, 2, 0, 1], C[1, 2, 1, 2], C[1, 2, 0, 2]],
                [C[0, 2, 0, 0], C[0, 2, 1, 1], C[0, 2, 2, 2], C[0, 2, 0, 1], C[0, 2, 1, 2], C[0, 2, 0, 2]],
            ]
        )


def undoVoigtMatrix(dim, C):
    """Put the material tangent moduli into matrix form.

    Parameters
    ----------
    dim
        Dimension of the model.
    C
        2D-Matrix to be put into 4D-matrix.

    Returns
    -------
    np.ndarray
        The strain in matrix notation for the corresponding dimension."""

    if dim == 2:
        return np.array(
            [
                [[[C[0, 0], C[0, 2]], [C[0, 2], C[0, 1]]], [[C[2, 0], C[2, 2]], [C[2, 2], C[2, 1]]]],
                [[[C[2, 0], C[2, 2]], [C[2, 2], C[2, 1]]], [[C[1, 0], C[1, 2]], [C[1, 2], C[1, 1]]]],
            ]
        )
    elif dim == 3:
        return np.array(
            [
                [
                    [[C[0, 0], C[0, 3], C[0, 5]], [C[0, 3], C[0, 1], C[0, 4]], [C[0, 5], C[0, 4], C[0, 2]]],
                    [[C[3, 0], C[3, 3], C[3, 5]], [C[3, 3], C[3, 1], C[3, 4]], [C[3, 5], C[3, 4], C[3, 2]]],
                    [[C[5, 0], C[5, 3], C[5, 5]], [C[5, 3], C[5, 1], C[5, 4]], [C[5, 5], C[5, 4], C[5, 2]]],
                ],
                [
                    [[C[3, 0], C[3, 3], C[3, 5]], [C[3, 3], C[3, 1], C[3, 4]], [C[3, 5], C[3, 4], C[3, 2]]],
                    [[C[1, 0], C[1, 3], C[1, 5]], [C[1, 3], C[1, 1], C[1, 4]], [C[1, 5], C[1, 4], C[1, 2]]],
                    [[C[4, 0], C[4, 3], C[4, 5]], [C[4, 3], C[4, 1], C[4, 4]], [C[4, 5], C[4, 4], C[4, 2]]],
                ],
                [
                    [[C[5, 0], C[5, 3], C[5, 5]], [C[5, 3], C[5, 1], C[5, 4]], [C[5, 5], C[5, 4], C[5, 2]]],
                    [[C[4, 0], C[4, 3], C[4, 5]], [C[4, 3], C[4, 1], C[4, 4]], [C[4, 5], C[4, 4], C[4, 2]]],
                    [[C[2, 0], C[2, 3], C[2, 5]], [C[2, 3], C[2, 1], C[2, 4]], [C[2, 5], C[2, 4], C[2, 2]]],
                ],
            ]
        )
