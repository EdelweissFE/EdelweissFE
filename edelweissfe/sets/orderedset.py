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
#  Paul Hofer Paul.Hofer@uibk.ac.at
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

from collections import UserDict
from collections.abc import Iterable
from types import MappingProxyType


class OrderedSet(UserDict):
    """An ordered set.

    Parameters
    ----------
    name
        A name for the object.
    items
        A list of items.
    """

    """ consider removing type checking! """

    def __init__(
        self,
        name: str,
        item_s,
    ):
        self.data = {}
        self.items = self.data.keys()
        self.keys = None
        self.values = None

        self.name = name

        # label is deprecated; use name instead!!
        self.label = name

        self.add(item_s)

    def checkObjectType(self, obj):
        """Checks if the object type is allowed in the OrderedSet"""

        return type(obj) in self.allowedObjectTypes

    def forceIter(self, item_s):
        """Return an iterator object for item_s even if item_s itself is not iterable"""

        if isinstance(item_s, Iterable):
            return iter(item_s)
        else:
            return iter([item_s])

    def add(self, item_s):
        """Add an item or an iterable of items to the OrderedSet"""

        # add items
        for item in self.forceIter(item_s):
            self.data.setdefault(item)

    def __setitem__(self, item, value):
        if self.checkObjectType(item):
            self.data[item] = None
        else:
            raise TypeError(f"You tried to add an item with wrong type: {item} of type {type(item)}")

    def __getitem__(self, key):
        return list(self.items)[key]

    # define &
    def __and__(self, other):
        if type(self) is type(other):
            return self.items & other.items
        else:
            raise TypeError("You can only compare OrderedSets with matching types")

    # define ^
    def __xor__(self, other):
        if type(self) is type(other):
            return self.items ^ other.items
        else:
            raise TypeError("You can only compare OrderedSets with matching types")

    # define |
    def __or__(self, other):
        if type(self) is type(other):
            return self.items | other.items
        else:
            raise TypeError("You can only compare OrderedSets with matching types")

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if type(self) is type(other):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __repr__(self):
        type_ = type(self)
        module = type_.__module__
        qualname = type_.__qualname__
        return (
            f'<{module}.{qualname} object at {hex(id(self))} with name "{self.name}">'
            # + f' "{self.name}"'
            # + "\n   ".join([""] + [item.__repr__() for item in self.data])
        )


class ImmutableOrderedSet(OrderedSet):
    """An immutable ordered set.

    Parameters
    ----------
    name
        A name for the object.
    items
        A list of items.
    """

    def __init__(
        self,
        name: str,
        item_s,
    ):
        filteredItems = dict()
        for item in self.forceIter(item_s):
            if self.checkObjectType(item):
                filteredItems.update({item: None})
            else:
                raise TypeError(f"You tried to add an item with wrong type: {item} of type {type(item)}")

        self.data = MappingProxyType(filteredItems)
        self.items = self.data.keys()
        self.keys = None
        self.values = None

        self.name = name
        self.name = name

    def add(self, item_s):
        raise TypeError(f"{type(self).__qualname__} items cannot be changed")

    def __setitem__(self, item, value):
        raise TypeError(f"{type(self).__qualname__} items cannot be changed")

    def __ior__(self, other):
        raise TypeError(f"{type(self).__qualname__} items cannot be changed")
