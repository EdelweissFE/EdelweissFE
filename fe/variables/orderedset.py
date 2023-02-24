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


class OrderedSet(list):
    """An ordered set.
    It has a label, and a list containing unique items.

    Parameters
    ----------
    label
        A label for the object.
    items
        A list of items.
    """

    allowedObjectTypes = []

    def __init__(
        self,
        label: str,
        items: list,
    ):
        self.label = label
        self.add(items)

    def checkObjectType(self, obj):
        """Checks if the object type is allowed for the respective OrderedSet

        Parameters
        ----------
        obj
            An arbitrary object to be checked
        """
        if type(obj) not in self.allowedObjectTypes:
            raise Exception("type {:} not allowed in {:}".format(type(obj), type(self)))
        return True

    def add(self, items: list):
        """Add a list of items to the OrderedSet.

        Parameters
        ----------
        items
            A list of items to be added.
        """
        super().extend([i for i in items if i not in self and self.checkObjectType(i)])

    def append(self, item):
        self.add([item])

    def extend(self, items: list):
        self.add(items)
