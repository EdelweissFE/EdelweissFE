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
# Created on Tue Jan  17 19:10:42 2017

# @author: Matthias Neuner
"""This section represents a classical plane solid materal section.
"""

documentation=
{"elementSets": "comma separated list of element sets for this section",
        "material": "the material to be assigned"}

import numpy as np
import gstools
from fe.utils.misc import stringDict


class Section:
    def __init__(self, name, options, materialName, thickness, modelInfo):

        self.materialName = materialName
        self.elSetNames = [e for l in options for e in l]
        self.thickness = thickness

    def assignSectionPropertiesToModel(self, modelInfo):

        elSets = [modelInfo["elementSets"][setName] for setName in self.elSetNames]
        material = modelInfo["materials"][self.materialName]

        for elSet in elSets:
            for el in elSet:
                elementThickness = self.thickness
                elProperties = np.array(
                    [
                        elementThickness,
                    ],
                    dtype=np.float,
                )
                el.setProperties(elProperties)
                el.initializeElement()
                el.setMaterial(material["name"], material["properties"])

        return modelInfo
