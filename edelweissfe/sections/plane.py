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

import numpy as np

from edelweissfe.sections.base.sectionbase import Section as SectionBase

"""This section represents a classical plane solid materal section.
"""

documentation = {
    "elementSets": "comma separated list of element sets for this section",
    "material": "the material to be assigned",
    "materialParameterFromField, index=[index of material parameter], value=[name of analytical field], type=[either 'setToValue' or 'scale']": "(optional) set or scale a material parameter using the value of the given analytical field; modify the field value using the optional keyword f(p,f)=[...] (p...value of parameter from material definition; f...value of analytical field)",
    "thickness": "the thickness to be assigned",
}


class Section(SectionBase):
    def __init__(self, name, options, materialName, model, **kwargs):
        super().__init__(name, options, materialName, model, **kwargs)
        try:
            self.thickness = kwargs["thickness"]
        except KeyError:
            raise KeyError(f"Thickness must be specified for section {name}")

    def assignSectionPropertiesToElement(self, element, **kwargs):
        material = kwargs.get("material", self.material)

        nSpatialDimensions = element.nSpatialDimensions
        if nSpatialDimensions != 2:
            raise Exception(f"Plane section is incompatible with {nSpatialDimensions}-dimensional finite elements.")

        thickness = self.thickness
        elProperties = np.array([thickness], dtype=float)

        element.setProperties(elProperties)
        element.initializeElement()
        # to make sure all elProviders work
        if not isinstance(material, dict):
            element.setMaterial(material)
        else:
            try:  # for Marmot
                element.setMaterial(material["name"], material["properties"])
            except TypeError:
                raise Exception("Material provider and element are not compatible!")
