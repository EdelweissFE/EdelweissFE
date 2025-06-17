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


class HyperplasticAdvancedMaterial(BaseHyperElasticMaterial):
    """Hyperplastic material with automatic differentiation for automatic calculation of the Kirchhoff stress and the tangent moduli.

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

        return 20

    def __init__(self, materialProperties: dict):
        self._materialProperties = materialProperties
        # elasticity parameters
        self._mu = materialProperties["mu"]
        self._K = materialProperties["K"]
        # plasticity parameters
        self.yieldStress = materialProperties["fy0"]
        self.HLin = materialProperties["HLin"]
        self.deltaYieldStress = materialProperties["dfy"]
        self.delta = materialProperties["delta"]
        self.setEnergyFunction(materialProperties["psi_e"])
        self._params = np.asarray(materialProperties.get("a", ""))  # extra optional parameters
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

        def getElasticEnergyAndStressAndTangent(Btr, dF, invFold, Bold):
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

            Returns
            -------
            double
                The energy density.
            np.ndarray
                The Kirchhoff stress.
            np.ndarray
                The elastic tangent dKirchhoff/dDeformationGradient."""

            J_ = np.sqrt(lin.det(Btr))
            I1_ = np.trace(Btr)
            invB = lin.inv(Btr)
            energy, dW_dI1, dW_dJ, d2W_dI1dJ = nd.second_partial_derivative(
                lambda I1, J: energyDensity(self._mu, self._K, I1, J, self._params), I1_, J_
            )
            d2W_dI12 = energyDensity(self._mu, self._K, nd.Dual2_64(I1_, 1.0, 0.0), J_, self._params).second_derivative
            d2W_dJ2 = energyDensity(self._mu, self._K, I1_, nd.Dual2_64(J_, 1.0, 0.0), self._params).second_derivative
            T = 2 * dW_dI1 * Btr + J_ * dW_dJ * Ieye
            dB_dF = np.einsum("ik,np,jp,ln->ijkl", Ieye, Bold, dF, invFold) + np.einsum(
                "io,on,jk,ln->ijkl", dF, Bold, Ieye, invFold
            )
            dJ_dF = 1 / 2 * J_ * np.einsum("nm,mnkl->kl", invB, dB_dF)
            dI1_dF = np.einsum("mn,mnkl->kl", Ieye, dB_dF)
            CTauF = (
                2 * (d2W_dI12 * outer2D(Btr, dI1_dF) + d2W_dI1dJ * outer2D(Btr, dJ_dF) + dW_dI1 * dB_dF)
                + (dW_dJ + J_ * d2W_dJ2) * outer2D(Ieye, dJ_dF)
                + J_ * d2W_dI1dJ * outer2D(Ieye, dI1_dF)
            )
            return (energy, T, CTauF)

        def stress(e):
            """Compute the Kirchhoff stress.

            Parameters
            ----------
            e
                The spatial Hencky strain.

            Returns
            -------
            np.ndarray
                The Kirchhoff stress."""

            J_ = np.exp(lin.trace(e))
            B, _ = tensorExp(2 * e)
            I1_ = np.trace(B)
            dW_dI1 = energyDensity(self._mu, self._K, nd.Dual64(I1_, 1.0), J_, self._params).first_derivative
            dW_dJ = energyDensity(self._mu, self._K, I1_, nd.Dual64(J_, 1.0), self._params).first_derivative
            return 2 * dW_dI1 * B + J_ * dW_dJ * Ieye

        def dKirchhoff_dE(e):
            """Compute the Kirchhoff stress.

            Parameters
            ----------
            e
                The spatial Hencky strain.

            Returns
            -------
            np.ndarray
                The derivative of the Kirchhoff stress w.r.t. the spatial Hencky strain."""

            J_ = np.exp(lin.trace(e))
            B, nExp = tensorExp(2 * e)
            I1_ = np.trace(B)
            dExpE = np.einsum("ijmn,mnkl->ijkl", dTensorExp_dA(2 * e, nExp), Ieye4D)
            _, dW_dI1, dW_dJ, d2W_dI1dJ = nd.second_partial_derivative(
                lambda I1, J: energyDensity(self._mu, self._K, I1, J, self._params), I1_, J_
            )
            d2W_dI12 = energyDensity(self._mu, self._K, nd.Dual2_64(I1_, 1.0, 0.0), J_, self._params).second_derivative
            d2W_dJ2 = energyDensity(self._mu, self._K, I1_, nd.Dual2_64(J_, 1.0, 0.0), self._params).second_derivative
            ExpP1 = np.einsum("ij,mn,mnkl->ijkl", B, Ieye, dExpE)
            ExpP2 = np.einsum("ij,mn,mnkl->ijkl", Ieye, Ieye, dExpE)
            IxI = outer2D(Ieye, Ieye)
            return (
                2 * d2W_dI12 * ExpP1
                + 2 * J_ * d2W_dI1dJ * outer2D(B, Ieye)
                + 4 * dW_dI1 * dExpE
                + J_ * dW_dJ * IxI
                + 2 * J_ * d2W_dI1dJ * ExpP2
                + J_**2 * d2W_dJ2 * IxI
            )

        self._getElasticResults = lambda Btr, dF, invFold, Bold: getElasticEnergyAndStressAndTangent(
            Btr, dF, invFold, Bold
        )
        self._kirchhoffStress = lambda e: stress(e)
        self._dKirchhoff_dE = lambda e: dKirchhoff_dE(e)

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
        # calculate elastic trial state
        (self._energy[0], kirchhoff, CTauF) = self._getElasticResults(Btrial, dF, invFold, Bold)
        # handle zero strain
        if np.all(np.abs(dF - Ieye) < 1e-14):
            dStress_dDeformationGradient[:] = CTauF
            return
        if self._f(deviatoricNorm(deviatoric(kirchhoff)), self._kappaOld) > 0.0:  # plastic
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
                # compute Kirchhoff stress
                kirchhoff = self._kirchhoffStress(e_e)
                devKirchhoff = deviatoric(kirchhoff)
                devNormKirchhoff = deviatoricNorm(devKirchhoff)
                n = devKirchhoff / devNormKirchhoff  # projection direction
                # compute the residuals
                Reps = doVoigtStrain(3, e_e - e_eTrial + np.sqrt(3 / 2) * dKappa * n)
                Rk = devNormKirchhoff - np.sqrt(2 / 3) * self._fy(self._kappaOld[0] + dKappa)
                # compute the residual derivations
                dT_de = doVoigtMatrix(3, np.einsum("ijkl,klmn->ijmn", IDev4D, self._dKirchhoff_dE(e_e))) @ Ihs
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
                dT_dE = doVoigtMatrix(3, self._dKirchhoff_dE(e)) @ Ihs
                return undoVoigtMatrix(3, dT_dE @ dX_dEtrial[:6, :6])

            # update to back projected state
            self._ee[:] = undoVoigtStrain(3, X[:6])
            dT_dE = elastoPlasticKirchhoffTangent(self._ee)
            # calculate back projected kirchhoff stress and convert to Voigt
            T = self._kirchhoffStress(self._ee)
            stress[:] = doVoigtStress(3, T)
            # compute elastoplastic consistent tangent
            dEe_dB = 1 / 2 * dTensorLog_dA(Btrial)
            dB_dF = np.einsum("ik,np,jp,ln->ijkl", Ieye, Bold, dF, invFold) + np.einsum(
                "io,on,jk,ln->ijkl", dF, Bold, Ieye, invFold
            )
            dEe_dF = np.einsum("ijmn,mnkl->ijkl", dEe_dB, dB_dF)
            dStress_dDeformationGradient[:] = np.einsum("ijop,opkl->ijkl", dT_dE, dEe_dF)
        else:  # elastic
            stress[:] = doVoigtStress(3, kirchhoff)
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
