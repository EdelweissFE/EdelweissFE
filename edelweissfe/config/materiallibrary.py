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

from edelweissfe.utils.misc import strCaseCmp


def getMaterialClass(materialName: str, provider: str = None) -> type:
    """Get the the requested material class.

    Parameters
    ----------
    materialName
        The name of the requested material.
    provider
        The name of the material provider.

    Returns
    -------
    type
        The material provider class type.
    """

    if provider is None:
        provider = "MarmotMaterial"

    if strCaseCmp(provider, "marmotmaterial"):

        return None

    if strCaseCmp(provider, "edelweiss"):
        if strCaseCmp(materialName, "linearelastic"):
            from edelweissfe.materials.linearelastic.linearelastic import (
                LinearElasticMaterial,
            )

            material = LinearElasticMaterial
        elif strCaseCmp(materialName, "vonmises"):
            from edelweissfe.materials.vonmises.vonmises import VonMisesMaterial

            material = VonMisesMaterial
        elif strCaseCmp(materialName, "neohookewa"):
            from edelweissfe.materials.neohooke.neohookepencegouformulationa import (
                NeoHookeanWaMaterial,
            )

            material = NeoHookeanWaMaterial
        elif strCaseCmp(materialName, "neohookewb"):
            from edelweissfe.materials.neohooke.neohookepencegouformulationb import (
                NeoHookeanWbMaterial,
            )

            material = NeoHookeanWbMaterial
        elif strCaseCmp(materialName, "neohookewc"):
            from edelweissfe.materials.neohooke.neohookepencegouformulationc import (
                NeoHookeanWcMaterial,
            )

            material = NeoHookeanWcMaterial
        elif strCaseCmp(materialName, "hyperelasticadvanced"):
            from edelweissfe.materials.hyperelasticadvanced.hyperelasticadvanced import (
                HyperelasticAdvancedMaterial,
            )

            material = HyperelasticAdvancedMaterial
        elif strCaseCmp(materialName, "hyperelasticadvancedi2extended"):
            from edelweissfe.materials.hyperelasticadvanced.hyperelasticadvancedi2extended import (
                HyperelasticAdvancedI2ExtendedMaterial,
            )

            material = HyperelasticAdvancedI2ExtendedMaterial
        elif strCaseCmp(materialName, "neohookewaplastic"):
            from edelweissfe.materials.neohookeplastic.neohookepencegouformulationaplastic import (
                NeoHookeanWaPlasticMaterial,
            )

            material = NeoHookeanWaPlasticMaterial
        elif strCaseCmp(materialName, "neohookewbplastic"):
            from edelweissfe.materials.neohookeplastic.neohookepencegouformulationbplastic import (
                NeoHookeanWbPlasticMaterial,
            )

            material = NeoHookeanWbPlasticMaterial
        elif strCaseCmp(materialName, "neohookewcplastic"):
            from edelweissfe.materials.neohookeplastic.neohookepencegouformulationcplastic import (
                NeoHookeanWcPlasticMaterial,
            )

            material = NeoHookeanWcPlasticMaterial
        elif strCaseCmp(materialName, "hyperplasticadvanced"):
            from edelweissfe.materials.hyperplasticadvanced.hyperplasticadvanced import (
                HyperplasticAdvancedMaterial,
            )

            material = HyperplasticAdvancedMaterial
        else:
            raise Exception("This material type doesn't exist (yet). Chosen material was: " + materialName)

        return material
