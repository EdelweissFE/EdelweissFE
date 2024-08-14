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

from edelweissfe.materials.base.basehypoelasticmaterial import BaseHypoElasticMaterial
from edelweissfe.utils.exceptions import CutbackRequest

# define variables and functions used in all material dimensions
IDev = np.array(
    [
        [2 / 3, -1 / 3, -1 / 3, 0, 0, 0],
        [-1 / 3, 2 / 3, -1 / 3, 0, 0, 0],
        [-1 / 3, -1 / 3, 2 / 3, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1],
    ]
)
IDevHalfShear = np.array(
    [
        [2 / 3, -1 / 3, -1 / 3, 0, 0, 0],
        [-1 / 3, 2 / 3, -1 / 3, 0, 0, 0],
        [-1 / 3, -1 / 3, 2 / 3, 0, 0, 0],
        [0, 0, 0, 0.5, 0, 0],
        [0, 0, 0, 0, 0.5, 0],
        [0, 0, 0, 0, 0, 0.5],
    ]
)
maxIter = 15  # max inner Newton cycles
tol = 1e-12  # tolerance


def deviatoric(stress):
    return IDev @ stress


def deviatoricNorm(devStress):
    return np.sqrt(np.sum(np.square(devStress[0:3])) + 2 * np.sum(np.square(devStress[3:6])))


class VonMisesMaterial(BaseHypoElasticMaterial):
    """Von Mises plastic material for 3D and 2D plane strain.

    Parameters
    ----------
    materialProperties
        The numpy array containing the material properties.

    Material properties
    -------------------
    E
        Elasticity module (Young's modulus).
    v
        Poisson's ratio.
    fy0
        Yield stress.
    HLin
        Linear plastic hardening parameter.
    dfy
        Multiplicator for nonlinear isotropic hardening.
    delta
        Exponent for nonlinear isotropic hardening.

    Hardening law
    -------------
    fy(kappa) = fy0 + HLin * kappa + dfy * (1 - exp(-delta * kappa))
    """

    def getNumberOfRequiredStateVars(self) -> int:
        """Returns number of needed material state Variables per integration point in the material.

        Returns
        -------
        int
            Number of needed material state Vars."""

        return 1

    def __init__(self, materialProperties: np.ndarray):
        # elasticity parameters
        self._E = materialProperties[0]  # set E
        self._v = materialProperties[1]  # set v
        # plasticity parameters
        self.yieldStress = materialProperties[2]
        self.HLin = materialProperties[3]
        self.deltaYieldStress = materialProperties[4]
        self.delta = materialProperties[5]
        # isotropic hardening
        self._fy = (
            lambda kappa_: self.yieldStress
            + self.HLin * kappa_
            + self.deltaYieldStress * (1.0 - np.exp(-self.delta * kappa_))
        )
        # yield function
        self._f = lambda devNorm, kappa_: devNorm - np.sqrt(2 / 3) * self._fy(kappa_)
        # derivative of fy dKappa
        self._dfy_ddKappa = lambda kappa_: self.HLin + self.deltaYieldStress * self.delta * np.exp(-self.delta * kappa_)
        self._G = self._E / (2 * (1.0 + self._v))

    def elasticityMatrix(self) -> np.ndarray:
        """Initalize a 3D material elasticity matrix.

        Returns
        -------
        np.ndarray
            The elasticity matrix."""

        # 3D elasticity matrix for Hexa8/Hexa20
        Ei = (
            self._E
            / ((1 + self._v) * (1 - 2 * self._v))
            * np.array(
                [
                    [(1 - self._v), self._v, self._v, 0, 0, 0],
                    [self._v, (1 - self._v), self._v, 0, 0, 0],
                    [self._v, self._v, (1 - self._v), 0, 0, 0],
                    [0, 0, 0, (1 - 2 * self._v) / 2, 0, 0],
                    [0, 0, 0, 0, (1 - 2 * self._v) / 2, 0],
                    [0, 0, 0, 0, 0, (1 - 2 * self._v) / 2],
                ]
            )
        )
        return Ei

    def assignCurrentStateVars(self, currentStateVars: np.ndarray):
        """Assign new current state vars.

        Parameters
        ----------
        currentStateVars
            Array containing the material state vars."""

        self.kappaOld = currentStateVars

    def computePlaneStress(
        self,
        stress: np.ndarray,
        dStressdStrain: np.ndarray,
        dStrain: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a plane stress material.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStressdStrain
            Matrix containing dStress/dStrain.
        dStrain
            Strain vector increment at time step t to t+dTime.
        time
            Array of step time and total time.
        dTime
            Current time step size."""

        raise Exception("Von Mises material with plane stress is currently not possible.")

    def computeStress(
        self,
        stress: np.ndarray,
        dStressdStrain: np.ndarray,
        dStrain: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a 3D material.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStressdStrain
            Matrix containing dStress/dStrain.
        dStrain
            Strain vector increment at time step t to t+dTime.
        time
            Array of step time and total time.
        dTime
            Current time step size."""

        Ei = self.elasticityMatrix()
        # handle zero strain increment
        if np.linalg.norm(dStrain) < 10**-14:
            dStressdStrain[:] = Ei
            return
        # elastic predictor
        trialStress = stress + dStrain @ Ei
        # deviatoric stress
        devStress = deviatoric(trialStress)
        # norm of deviatoric stress
        devNorm_ = deviatoricNorm(devStress)
        if self._f(devNorm_, self.kappaOld[0]) > 0.0:  # plastic step
            R = (
                lambda deltaKappa: devNorm_
                - np.sqrt(6) * self._G * deltaKappa
                - np.sqrt(2 / 3) * self._fy(self.kappaOld[0] + deltaKappa)
            )
            n = devStress / devNorm_  # projection direction
            counter = 0
            dKappa = 0
            while np.abs(R(dKappa)) > tol:
                if counter == maxIter:
                    raise CutbackRequest("Von Mises Newton failed.", 0.5)
                # compute derivative of R with dKappa
                dR_ddKappa = -np.sqrt(6) * self._G - np.sqrt(2 / 3) * self._dfy_ddKappa(self.kappaOld[0] + dKappa)
                # update dKappa and iteration counter
                dKappa -= R(dKappa) / dR_ddKappa
                counter += 1
            dLambda = np.sqrt(3 / 2) * dKappa
            # update kappa
            self.kappaOld[0] += dKappa
            stress[:] = trialStress - 2.0 * self._G * dLambda * n
            dStressdStrain[:] = (
                Ei
                - 2.0
                * self._G
                * (
                    1.0 / (1.0 + self._dfy_ddKappa(self.kappaOld[0]) / (3.0 * self._G))
                    - 2.0 * self._G * dLambda / devNorm_
                )
                * np.outer(n, n)
                - 4.0 * self._G**2 * dLambda / devNorm_ * IDevHalfShear
            )
        else:  # elastic step
            stress[:] = trialStress
            dStressdStrain[:] = Ei

    def computeUniaxialStress(
        self,
        stress: np.ndarray,
        dStressdStrain: np.ndarray,
        dStrain: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a uniaxial stress state.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStressdStrain
            Matrix containing dStress/dStrain.
        dStrain
            Strain vector increment at time step t to t+dTime.
        time
            Array of step time and total time.
        dTime
            Current time step size."""

        raise Exception("Computing uniaxial stress is not possible with this material.")

    def getResult(self, result: str) -> float:
        """Get the result, as a persistent view which is continiously
        updated by the material.

        Parameters
        ----------
        result
            The name of the result.

        Returns
        -------
        float
            The result.
        """

        if result == "kappa":
            return self.kappaOld
        else:
            raise Exception("This result doesn't exist for the current material.")
