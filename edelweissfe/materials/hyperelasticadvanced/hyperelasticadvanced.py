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

import num_dual as nd
import numpy as np
import numpy.linalg as lin

from edelweissfe.materials.base.basehyperelasticmaterial import BaseHyperElasticMaterial
from edelweissfe.utils.voigtnotation import doVoigtStress

# define identity matrix
Ieye = np.eye(3)


def outer2D(A, B):
    """Computes the outer product of two matrices, the result is a 4D matrix.

    Parameters
    ----------
    A
        First matrix to compute the outer product with.
    B
        Second matrix to compute the outer product with.

    Returns
    -------
    np.ndarray
        The 4D outer product."""

    return np.einsum("ij,kl->ijkl", A, B)


class HyperelasticAdvancedMaterial(BaseHyperElasticMaterial):
    """Hyperelastic material with automatic differentiation for automatic calculation of the Kirchhoff stress and the tangent moduli.

    Parameters
    ----------
    materialProperties
        The dictionary containing the material properties for the requested material."""

    @property
    def materialProperties(self) -> dict:
        """The properties the material has."""

        return self._materialProperties

    def getNumberOfRequiredStateVars(self) -> int:
        """Returns number of needed material state Variables per integration point in the material.

        Returns
        -------
        int
            Number of needed material state Vars."""

        return 1

    def __init__(self, materialProperties: dict):
        self._materialProperties = materialProperties
        # elastic parameters
        self._mu = materialProperties["mu"]
        self._K = materialProperties["K"]
        self.setEnergyFunction(materialProperties["psi_e"])
        self._params = np.asarray(materialProperties.get("a", ""))  # extra optional parameters

    def setEnergyFunction(self, materialEnergy: str):
        """Sets the energy density function for the custom material.

        Parameters
        ----------
        materialEnergy
            Energy density function for the material."""

        self._materialEnergy = materialEnergy

        if materialEnergy is None:
            raise Exception("This material requires an energy density function set with W='fct(C)'.")

        def energyDensity(mu, K, I1, J, a):
            """The energy density function.

            Parameters
            ----------
            mu
                The shear modulus.
            K
                The bulk modulus.
            I1
                The first invariant of a Cauchy-Green tensor.
            J
                The third invariant of a Cauchy-Green tensor.
            a
                User-given custom constant parameters.

            Returns
            -------
            double
                The energy density."""

            try:
                return eval(materialEnergy)
            except NameError:
                raise Exception(
                    "This material only allows parameters mu, K, the first invariant of C (I1) and J = det(F) as well as your own parameters a[i]."
                )

        def getElasticEnergyAndStressAndTangent(F):  # energy, kirchhoff stress and tangent
            """Compute the energy, the Kirchhoff stress and the tangent.

            Parameters
            ----------
            F
                The deformation gradient.

            Returns
            -------
            double
                The energy density.
            np.ndarray
                The Kirchhoff stress.
            np.ndarray
                The tangent dKirchhoff/dDeformationGradient."""

            B = F @ F.T
            J_ = np.sqrt(lin.det(B))
            I1_ = np.trace(B)
            invF = lin.inv(F)
            energy, dW_dI1, dW_dJ, d2W_dI1dJ = nd.second_partial_derivative(
                lambda I1, J: energyDensity(self._mu, self._K, I1, J, self._params), I1_, J_
            )
            d2W_dI12 = energyDensity(self._mu, self._K, nd.Dual2_64(I1_, 1.0, 0.0), J_, self._params).second_derivative
            d2W_dJ2 = energyDensity(self._mu, self._K, I1_, nd.Dual2_64(J_, 1.0, 0.0), self._params).second_derivative
            T = 2 * dW_dI1 * B + J_ * dW_dJ * Ieye
            dB_dF = np.einsum("ik,jl->ijkl", Ieye, F) + np.einsum("il,jk->ijkl", F, Ieye)
            CTauF = (
                2 * (2 * d2W_dI12 * outer2D(B, F) + J_ * d2W_dI1dJ * outer2D(B, invF.T) + dW_dI1 * dB_dF)
                + (J_ * dW_dJ + J_**2 * d2W_dJ2) * outer2D(Ieye, invF.T)
                + 2 * J_ * d2W_dI1dJ * outer2D(Ieye, F)
            )
            return (energy, T, CTauF)

        self._getElasticResults = lambda F: getElasticEnergyAndStressAndTangent(F)

    def assignCurrentStateVars(self, currentStateVars: np.ndarray):
        """Assign new current state vars.

        Parameters
        ----------
        currentStateVars
            Array containing the material state vars."""

        self._energy = currentStateVars

    def computePlaneKirchhoff(
        self,
        stress: np.ndarray,
        dStress_dDeformationGradient: np.ndarray,
        deformationGradient: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a 2D material with plane stress.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStress_dDeformationGradient
            Matrix containing dStress/dStrain.
        deformationGradient
            The deformation gradient at time step t.
        time
            Array of step time and total time.
        dTime
            Current time step size."""

        super().computePlaneKirchhoff(stress, dStress_dDeformationGradient, deformationGradient, time, dTime)

    def computeKirchhoff(
        self,
        stress: np.ndarray,
        dStress_dDeformationGradient: np.ndarray,
        deformationGradient: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a 3D material/2D material with plane strain.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStress_dDeformationGradient
            Matrix containing dStress/dStrain.
        deformationGradient
            The deformation gradient at time step t.
        time
            Array of step time and total time.
        dTime
            Current time step size."""

        (self._energy[0], T, CTauF) = self._getElasticResults(deformationGradient)
        stress[:] = doVoigtStress(3, T)
        dStress_dDeformationGradient[:] = CTauF

    def computeUniaxialKirchhoff(
        self,
        stress: np.ndarray,
        dStress_dDeformationGradient: np.ndarray,
        deformationGradient: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a uniaxial stress state.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStress_dDeformationGradient
            Matrix containing dStress/dStrain.
        deformationGradient
            The deformation gradient at time step t.
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

        if result == "energy":
            return self._energy
        else:
            raise Exception("This result doesn't exist for the current material.")
