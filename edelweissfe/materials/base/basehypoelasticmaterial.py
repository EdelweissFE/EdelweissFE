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

from abc import ABC, abstractmethod

import numpy as np

from edelweissfe.utils.exceptions import CutbackRequest


class BaseHypoElasticMaterial(ABC):
    """Base material class for a hypo elastic material.

    Parameters
    ----------
    materialProperties
        The numpy array containing the material properties for the requested material."""

    @property
    def materialProperties(self) -> np.ndarray:
        """The properties the material has."""

    @abstractmethod
    def getNumberOfRequiredStateVars(self) -> int:
        """Returns number of needed material state Variables per integration point in the material.

        Returns
        -------
        int
            Number of needed material state Vars."""

    @abstractmethod
    def __init__(self, materialProperties: np.ndarray):
        """Initialize."""

    @abstractmethod
    def assignCurrentStateVars(self, currentStateVars: np.ndarray):
        """Assign new current state vars.

        Parameters
        ----------
        currentStateVars
            Array containing the material state vars."""

    @abstractmethod
    def computePlaneStress(
        self,
        stress: np.ndarray,
        dStress_dStrain: np.ndarray,
        dStrain: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a 2D material with plane stress.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStress_dStrain
            Matrix containing dStress/dStrain.
        dStrain
            Strain vector increment at time step t to t+dTime.
        time
            Array of step time and total time.
        dTime
            Current time step size."""

        dStress_dStrain3D = np.zeros([6, 6])
        dStrain[2] = -dStrain[0] - dStrain[1]
        counter = 1
        while True:
            stress3D = np.zeros([6])
            self.computeStress(stress3D, dStress_dStrain3D, dStrain, time, dTime)
            residual = abs(stress3D[2])  # make S33 the residual
            if residual < 1e-11 or (counter > 7 and residual < 1e-8):
                break
            revTangent = 1.0 / dStress_dStrain3D[2, 2]
            if np.isnan(revTangent) or np.abs(revTangent) > 1e10:
                revTangent = 1e10 * np.sign(revTangent)
            dStrain[2] -= revTangent * stress3D[2]
            counter += 1
            if counter > 13:
                raise CutbackRequest("Plane Stress Newton failed.", 0.5)
        stress[:] += stress3D
        C = dStress_dStrain3D
        dStress_dStrain[:] = np.array(
            [
                [
                    C[0, 0] - C[2, 0] * C[0, 2] / C[2, 2],
                    C[0, 1] - C[2, 1] * C[0, 2] / C[2, 2],
                    C[0, 3] - C[2, 3] * C[0, 2] / C[2, 2],
                ],
                [
                    C[1, 0] - C[2, 0] * C[1, 2] / C[2, 2],
                    C[1, 1] - C[2, 1] * C[1, 2] / C[2, 2],
                    C[1, 3] - C[2, 3] * C[1, 2] / C[2, 2],
                ],
                [
                    C[3, 0] - C[2, 0] * C[3, 2] / C[2, 2],
                    C[3, 1] - C[2, 1] * C[3, 2] / C[2, 2],
                    C[3, 3] - C[2, 3] * C[3, 2] / C[2, 2],
                ],
            ]
        )

    @abstractmethod
    def computeStress(
        self,
        stress: np.ndarray,
        dStress_dStrain: np.ndarray,
        dStrain: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a 3D material/2D material with plane strain.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStress_dStrain
            Matrix containing dStress/dStrain.
        dStrain
            Strain vector increment at time step t to t+dTime.
        time
            Array of step time and total time.
        dTime
            Current time step size."""

    @abstractmethod
    def computeUniaxialStress(
        self,
        stress: np.ndarray,
        dStress_dStrain: np.ndarray,
        dStrain: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a uniaxial stress state.

        Parameters
        ----------
        stress
            Vector containing the stresses.
        dStress_dStrain
            Matrix containing dStress/dStrain.
        dStrain
            Strain vector increment at time step t to t+dTime.
        time
            Array of step time and total time.
        dTime
            Current time step size."""

    @abstractmethod
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
