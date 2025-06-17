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

from edelweissfe.materials.base.basehyperelasticmaterial import BaseHyperElasticMaterial
from edelweissfe.utils.exceptions import CutbackRequest
from edelweissfe.utils.plasticityutils import (
    deviatoric,
    deviatoricNorm,
    dTensorExp_dA,
    dTensorLog_dA,
    tensorExp,
    tensorLogarithmEig,
)
from edelweissfe.utils.voigtnotation import (
    doVoigtMatrix,
    doVoigtStrain,
    doVoigtStress,
    undoVoigtMatrix,
    undoVoigtStrain,
)

# define variables and functions
maxIter = 20  # max inner Newton cycles
tol = 1e-12  # tolerance

# deviatoric multiplicator
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
IDev4D = undoVoigtMatrix(3, IDev)
Ihs = np.diag([1, 1, 1, 0.5, 0.5, 0.5])
Ids = np.diag([1, 1, 1, 2, 2, 2])
Ieye = np.eye(3)
Ieye4D = undoVoigtMatrix(3, np.eye(6))
dR_dETrial = np.vstack([-np.eye(6), [np.zeros([6])]])


# material-specific functions
def energyAndStressAndElasticTangent(Btr, dF, invFold, Bold, K, mu):
    """Compute the energy, the Kirchhoff stress and the elastic tangent.

    Parameters
    ----------
    Btr
        The trial left Cauchy-Green tensor.
    dF
        The change of the deformation gradient since the last step.
    invFold
        The inverse of the deformation gradient from the last step.
    Bold
        The left Cauchy-Green tensor from the last step.
    K
        The bulk modulus.
    mu
        The shear modulus.

    Returns
    -------
    double
        The energy density.
    np.ndarray
        The Kirchhoff stress.
    np.ndarray
        The elastic tangent dKirchhoff/dDeformationGradient."""

    J = np.sqrt(lin.det(Btr))
    I1 = np.trace(Btr)
    invB = lin.inv(Btr)
    muBar = mu * J ** (2 / 3 - K / mu)
    lambdaBar = (K / mu - 2 / 3) * muBar / 2
    T = mu * Btr - muBar * np.eye(3)
    dB_dF = np.einsum("ik,np,jp,ln->ijkl", Ieye, Bold, dF, invFold) + np.einsum(
        "io,on,jk,ln->ijkl", dF, Bold, Ieye, invFold
    )
    CTauF = mu * dB_dF + lambdaBar * np.einsum("ij,mn,mnkl->ijkl", Ieye, invB.T, dB_dF)
    energy = mu / 2 * (I1 - 3) + 3 * mu**2 / (2 * K - 2 * mu) * (J ** (2 / 3 - K / mu) - 1)
    return (energy, T, CTauF)


def kirchhoffStress(e, K, mu):
    """Compute the derivative of the Kirchhoff stress w.r.t. the spatial logarithmic Hencky strain.

    Parameters
    ----------
    e
        The spatial logarithmic Hencky strain.
    K
        The bulk modulus.
    mu
        The shear modulus.

    Returns
    -------
    np.ndarray
        The Kirchhoff stress."""

    B, _ = tensorExp(2 * e)
    muBar = np.exp((2 / 3 - K / mu) * lin.trace(e))
    T = mu * (B - muBar * Ieye)
    return T


def dKirchhoff_dE(e, K, mu):
    """Compute the derivative of the Kirchhoff stress w.r.t. the spatial logarithmic Hencky strain.

    Parameters
    ----------
    e
        The spatial logarithmic Hencky strain.
    K
        The bulk modulus.
    mu
        The shear modulus.

    Returns
    -------
    np.ndarray
        The derivative of the Kirchhoff stress w.r.t. the spatial Hencky strain."""

    B, n = tensorExp(2 * e)
    J = np.sqrt(lin.det(B))
    dExpE = np.einsum("ijmn,mnkl->ijkl", dTensorExp_dA(2 * e, n), Ieye4D)
    return 2 * mu * dExpE + (K - 2 / 3 * mu) * J ** (2 / 3 - K / mu) * np.einsum("ij,kl->ijkl", Ieye, Ieye)


class NeoHookeanWcPlasticMaterial(BaseHyperElasticMaterial):
    """Neo-Hookean Wb material according to [1] with the following energy density function.

    W_c (I1, J) = mu/2 * (I1 - 3) + 3mu²/(3K - 2mu) * (J^(2/3 - K/mu) - 1).

    [1] Pence, T. J., & Gou, K. (2015). On compressible versions of the incompressible neo-Hookean material.
        Mathematics and Mechanics of Solids, 20(2), 157–182. https://doi.org/10.1177/1081286514544258

    Parameters
    ----------
    materialProperties
        The numpy array containing the material properties for the requested material."""

    @property
    def materialProperties(self) -> np.ndarray:
        """The properties the material has."""

        return self._materialProperties

    def getNumberOfRequiredStateVars(self) -> int:
        """Returns number of needed material state Variables per integration point in the material.

        Returns
        -------
        int
            Number of needed material state Vars."""

        return 20

    def __init__(self, materialProperties: np.ndarray):
        self._materialProperties = materialProperties
        # elasticity parameters
        self._mu = materialProperties[0]
        self._K = materialProperties[1]
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

    def assignCurrentStateVars(self, currentStateVars: np.ndarray):
        """Assign new current state vars.

        Parameters
        ----------
        currentStateVars
            Array containing the material state vars."""

        self._energy = np.reshape(currentStateVars[0:1], (1))
        self._kappaOld = np.reshape(currentStateVars[1:2], (1))
        self._Fold = np.reshape(currentStateVars[2:11], (3, 3))
        self._ee = np.reshape(currentStateVars[11:], (3, 3))
        if np.all(self._Fold == 0):
            self._Fold[:] = np.eye(3)

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
        invFold = lin.inv(self._Fold)
        dF = F @ invFold
        # update to new values
        Bold, _ = tensorExp(2 * self._ee)
        Btrial = dF @ Bold @ dF.T
        e_eTrial = 1 / 2 * tensorLogarithmEig(Btrial)
        self._energy[0], trialKirchhoff, CTauF = energyAndStressAndElasticTangent(
            Btrial, dF, invFold, Bold, self._K, self._mu
        )
        # handle zero strain
        if np.all(np.abs(dF - Ieye) < 1e-14):
            dStress_dDeformationGradient[:] = CTauF
            return
        # calculate deviatoric norm of Kirchhoff stress
        devNormKirchhoff = deviatoricNorm(deviatoric(trialKirchhoff))
        if self._f(devNormKirchhoff, self._kappaOld[0]) > 0.0:  # plastic step
            # small strain algorithm to solve for kappa and epsilon

            def residuals(X):
                """Compute the residual and its differential for the plastic solving process.

                Parameters
                ----------
                X
                    The current solution vector.

                Returns
                -------
                np.ndarray
                    The residual vector.
                np.ndarray
                    The derivative of the residual vector w.r.t. X."""

                e_e = undoVoigtStrain(3, X[:6])
                dKappa = X[6]
                # compute Kirchhoff stress and Tensor exponential row number
                kirchhoff = kirchhoffStress(e_e, self._K, self._mu)
                devKirchhoff = deviatoric(kirchhoff)
                devNormKirchhoff = deviatoricNorm(devKirchhoff)
                n = devKirchhoff / devNormKirchhoff  # projection direction
                # compute the residuals
                Reps = doVoigtStrain(3, e_e - e_eTrial + np.sqrt(3 / 2) * dKappa * n)
                Rk = devNormKirchhoff - np.sqrt(2 / 3) * self._fy(self._kappaOld[0] + dKappa)
                # compute the residual derivations
                dT_de = (
                    doVoigtMatrix(3, np.einsum("ijkl,klmn->ijmn", IDev4D, dKirchhoff_dE(e_e, self._K, self._mu))) @ Ihs
                )
                nv = doVoigtStress(3, n)  # make Voigt
                dReps_dEe = np.eye(6) + np.sqrt(3 / 2) * dKappa / devNormKirchhoff * (
                    dT_de - Ids @ np.outer(nv, nv) @ dT_de
                )
                dReps_dKappa = Ids @ np.vstack(np.sqrt(3 / 2) * nv)
                dRk_dKappa = -np.sqrt(2 / 3) * self._dfy_ddKappa(self._kappaOld[0] + dKappa)
                # combine derivations into one matrix
                dR_dX = np.vstack([np.hstack([dReps_dEe, dReps_dKappa]), np.hstack([nv @ dT_de, [dRk_dKappa]])])
                return (np.hstack([Reps, [Rk]]), dR_dX)

            # start values
            X = np.hstack([doVoigtStrain(3, e_eTrial), 0])
            counter = 0
            R, dR_dX = residuals(X)
            relTol = tol * self._fy(self._kappaOld[0])
            while lin.norm(R) > relTol:  # check if residuals fulfill tolerances
                if counter == maxIter:
                    raise CutbackRequest("J2 plasticity solving failed.", 0.5)
                # compute residual and its derivation
                R, dR_dX = residuals(X)
                X -= lin.solve(dR_dX, R)
                counter += 1
            # update kappa
            self._kappaOld[:] += X[6]

            # calculation of the elastoplastic consistent tangent
            def elastoPlasticKirchhoffTangent(e):
                """Compute the elastoplastic Kirchhoff tangent dKirchhoff/dEtrial.

                Parameters
                ----------
                e
                    The elastic spatial Hencky strain.

                Returns
                -------
                np.ndarray
                    The elastoplastic tangent dKirchhoff/dEtrial."""

                dX_dEtrial = -lin.solve(dR_dX, dR_dETrial)
                dT_dE = doVoigtMatrix(3, dKirchhoff_dE(e, self._K, self._mu)) @ Ihs
                return undoVoigtMatrix(3, dT_dE @ dX_dEtrial[:6, :6])

            # update to back projected state
            self._ee[:] = undoVoigtStrain(3, X[:6])
            dT_dE = elastoPlasticKirchhoffTangent(self._ee)
            # calculate back projected kirchhoff stress and convert to Voigt
            T = kirchhoffStress(self._ee, self._K, self._mu)
            stress[:] = doVoigtStress(3, T)
            # compute elastoplastic consistent tangent
            dEe_dB = 1 / 2 * dTensorLog_dA(Btrial)
            dB_dF = np.einsum("ik,np,jp,ln->ijkl", Ieye, Bold, dF, invFold) + np.einsum(
                "io,on,jk,ln->ijkl", dF, Bold, Ieye, invFold
            )
            dEe_dF = np.einsum("ijmn,mnkl->ijkl", dEe_dB, dB_dF)
            dStress_dDeformationGradient[:] = np.einsum("ijop,opkl->ijkl", dT_dE, dEe_dF)
        else:  # elastic step
            stress[:] = doVoigtStress(3, trialKirchhoff)
            dStress_dDeformationGradient[:] = CTauF
            self._ee[:] = e_eTrial
        self._Fold[:] = F

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
        elif result == "kappa":
            return self._kappaOld
        else:
            raise Exception("This result doesn't exist for the current material.")
