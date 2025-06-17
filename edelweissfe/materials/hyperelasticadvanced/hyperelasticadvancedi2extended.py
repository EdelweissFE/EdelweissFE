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

from edelweissfe.materials.base.basehyperelasticmaterial import BaseHyperElasticMaterial
from edelweissfe.utils.voigtnotation import (
    doVoigtMatrix,
    doVoigtStress,
    undoVoigtMatrix,
    undoVoigtStress,
)


def putCInVoigt(C):
    """Put the right Cauchy-Green tensor C in length 9 Voigt notation.

    Parameters
    ----------
    C
        The right Cauchy-Green tensor."""

    return np.array([C[0, 0], C[1, 1], C[2, 2], C[0, 1], C[1, 2], C[0, 2], C[1, 0], C[2, 1], C[2, 0]])


def I1(C):  # trace of C in Voigt
    """Compute the first invariant with a vectorized C.

    Parameters
    ----------
    C
        The right Cauchy-Green tensor in length 9 Voigt notation.

    Returns
    -------
    double
        The first invariant of the right Cauchy-Green tensor."""

    return C[0] + C[1] + C[2]


def I2(C):  # second invariant in Voigt
    """Compute the second invariant with a vectorized C.

    Parameters
    ----------
    C
        The right Cauchy-Green tensor in length 9 Voigt notation.

    Returns
    -------
    double
        The second invariant of the right Cauchy-Green tensor."""

    Cm = np.array([[C[0], C[3], C[5]], [C[6], C[1], C[4]], [C[8], C[7], C[2]]])
    CC = Cm @ Cm
    CCv = putCInVoigt(CC)
    return 1 / 2 * (I1(C) ** 2 - I1(CCv))


def J(C):  # det of F in Voigt
    """Compute the determinant of the deformation gradient with a vectorized C.

    Parameters
    ----------
    C
        The right Cauchy-Green tensor in length 9 Voigt notation.

    Returns
    -------
    double
        The determinant of the deformation gradient."""

    return np.sqrt(
        C[0] * C[1] * C[2]
        + C[3] * C[4] * C[5]
        + C[6] * C[7] * C[8]
        - C[1] * C[5] * C[8]
        - C[0] * C[4] * C[7]
        - C[2] * C[3] * C[6]
    )


class HyperelasticAdvancedI2ExtendedMaterial(BaseHyperElasticMaterial):
    """Hyperelastic material with automatic differentiation for automatic calculation of the Kirchhoff stress and the tangent moduli.

    Parameters
    ----------
    materialProperties
        The numpy array containing the material properties for the requested material."""

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
        if "C" not in materialEnergy:  # use num_dual
            import num_dual as nd
            import numpy as np

            self._autodiffisnum_dual = True
        else:  # use autograd
            import autograd as ag
            import autograd.numpy as na

            self._autodiffisnum_dual = False

        def energyDensity(C, mu, K, I1, I2, J, a, np):
            """The energy density function.

            Parameters
            ----------
            C
                The right Cauchy-Green tensor.
            mu
                The shear modulus.
            K
                The bulk modulus.
            I1
                The first invariant of a Cauchy-Green tensor.
            I2
                The second invariant of a Cauchy-Green tensor.
            J
                The third invariant of a Cauchy-Green tensor.
            a
                User-given custom constant parameters.
            np
                The local package implementation of numpy.

            Returns
            -------
            double
                The energy density."""

            try:
                return eval(materialEnergy)
            except NameError:
                raise Exception(
                    "This material only allows parameters mu, K, the first 2 invariants of C (I1, I2), C itself and J = det(F) as well as your own parameters a[i]."
                )

        def getEnergyAndStressAndTangent(C, energyFct):
            """Compute the energy, the Kirchhoff stress and the tangent.

            Parameters
            ----------
            C
                The right Cauchy-Green tensor.
            energyFct
                The user-given function for the energy density.

            Returns
            -------
            double
                The energy density.
            np.ndarray
                The Kirchhoff stress.
            np.ndarray
                The tangent dKirchhoff/dDeformationGradient."""

            if self._autodiffisnum_dual:
                Cv = putCInVoigt(C)
                (energy, stressF, CSEF) = nd.hessian(
                    lambda Cv: energyFct(Cv, self._mu, self._K, I1(Cv), I2(Cv), J(Cv), self._params, np), Cv
                )
                CSEF = np.array(CSEF)
                CSE = 2 * CSEF[:6, :6]
                CSE[:3, :3] = 2 * CSE[:3, :3]
                CSE[3:, :3] += 2 * CSEF[6:, :3]
                CSE[:3, 3:] += 2 * CSEF[:3, 6:]
                CSE[3:, 3:] = CSEF[3:6, 3:6] + CSEF[6:, 3:6] + CSEF[3:6, 6:] + CSEF[6:, 6:]
                return (energy, 2 * np.array(stressF[:6]), CSE)
            else:

                def W(C):
                    return energyDensity(
                        C,
                        self._mu,
                        self._K,
                        na.trace(C),
                        1 / 2 * ((na.trace(C)) ** 2 - na.trace(C @ C)),
                        na.sqrt(na.linalg.det(C)),
                        self._params,
                        na,
                    )

                def S(C):
                    return ag.grad(W)(C) + ag.grad(lambda C: W(C.T))(C)

                def CSE(C):
                    return ag.jacobian(S)(C) + ag.jacobian(lambda C: S(C.T))(C)

                return (W(C), doVoigtStress(3, S(C)), doVoigtMatrix(3, CSE(C)))

        self._getHessian = lambda C: getEnergyAndStressAndTangent(C, energyDensity)

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

        F = deformationGradient
        (self._energy[0], S, CSE) = self._getHessian(F.T @ F)
        S = undoVoigtStress(3, S)
        stress[:] = doVoigtStress(3, F @ S @ F.T)
        dC_dF = np.einsum("lo,kp->opkl", np.eye(3), F) + np.einsum("ko,lp->opkl", F, np.eye(3))
        CTauF1 = np.einsum("ik,ln,jn->ijkl", np.eye(3), S, F)
        CTauF2 = np.einsum("im,mnop,opkl,jn->ijkl", F, 1 / 2 * undoVoigtMatrix(3, CSE), dC_dF, F)
        CTauF3 = np.einsum("im,ml,jk->ijkl", F, S, np.eye(3))
        dStress_dDeformationGradient[:] = CTauF1 + CTauF2 + CTauF3

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
