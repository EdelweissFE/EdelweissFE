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

from abc import ABC, abstractmethod

import numpy as np


class BaseHypoElasticMaterial(ABC):
    """Base material class for a hypo elastic material.

    Parameters
    ----------
    materialProperties
        The numpy array containing the material properties for the requested material."""

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
        dStressdStrain: np.ndarray,
        dStrain: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a 2D material with plane stress.

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

    @abstractmethod
    def computeStress(
        self,
        stress: np.ndarray,
        dStressdStrain: np.ndarray,
        dStrain: np.ndarray,
        time: float,
        dTime: float,
    ):
        """Computes the stresses for a 3D material/2D material with plane strain.

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

    @abstractmethod
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