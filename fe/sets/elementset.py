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
#  Alexander Dummer alexander.dummer@uibk.ac.at
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

# @author: Alexander Dummer

from fe.sets.orderedset import OrderedSet
from fe.elements.marmotelement.element import MarmotElementWrapper
from fe.elements.marmotsingleqpelement.element import MarmotMaterialWrappingElement


class ElementSet(OrderedSet):
    """A basic element set.
    It has a label, and a list containing unique elements.

    Parameters
    ----------
    label
        The unique label for this element set.
    elements
        A list of elements.
    """

    allowedObjectTypes = [MarmotElementWrapper, MarmotMaterialWrappingElement]

    def __init__(
        self,
        label: str,
        elements: list,
    ):
        super().__init__(label, elements)
