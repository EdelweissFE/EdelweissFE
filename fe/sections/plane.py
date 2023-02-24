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
#  Paul Hofer paul.hofer@uibk.ac.at
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
#
# @author: Matthias Neuner, Paul Hofer

"""This section represents a classical plane solid materal section.
"""

documentation = {
    "elementSets": "comma separated list of element sets for this section",
    "material": "the material to be assigned",
    "materialParameterFromField, index=[index of material parameter], value=[name of analytical field], type=[either 'setToValue' or 'scale']": "(optional) set or scale a material parameter using the value of the given analytical field; modify the field value using the optional keyword f(p,f)=[...] (p...value of parameter from material definition; f...value of analytical field)",
    "thickness": "the thickness to be assigned",
}

import numpy as np

from fe.sections.base.sectionbase import Section as SectionBase


class Section(SectionBase):
    def __init__(self, name, options, materialName, thickness, model):
        super().__init__(name, options, materialName, thickness, model)
        self.thickness = thickness
        if not thickness > 0:
            raise ValueError("{}: Thickness must be greater than zero".format(name))

    def assignSectionPropertiesToModel(self, model):
        elSets = [model["elementSets"][setName] for setName in self.elSetNames]
        material = model["materials"][self.materialName]

        elProperties = np.array(
            [
                self.thickness,
            ],
            dtype=float,
        )

        if any(self.materialParameterFromFieldDefs):
            for elSet in elSets:
                for el in elSet:
                    materialProperties = super().propertiesFromField(el, material, model)

                    el.setProperties(elProperties)
                    el.initializeElement()
                    el.setMaterial(material["name"], materialProperties)

        else:
            for elSet in elSets:
                for el in elSet:
                    el.setProperties(elProperties)
                    el.initializeElement()
                    el.setMaterial(material["name"], material["properties"])

        if self.writeMaterialPropertiesToFile:
            self.exportMaterialPropertiesToFile(elSets)

        return model
