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
"""
EdelweissFE currently supports finite element implementations provided by the Marmot library.
In future, elements by other providers or elements directly implemented in EdelweissFE may be added here.

.. code-block:: edelweiss
    :caption: Example:

    *element, type=C3D8, provider=marmot
        ** el_label, node1, node2, node3, node4, ...
        1000,        1,     2,     3,     4,     ...
"""

from edelweissfe.elements.displacementelement.element import (  # noqa: F401
    DisplacementElement,
)
from edelweissfe.elements.displacementtlelement.element import (  # noqa: F401
    DisplacementTLElement,
)
from edelweissfe.elements.library import elLibrary
from edelweissfe.utils.misc import strCaseCmp


def getElementClass(elType: str, provider: str = None) -> type:
    """Get the class type of the requested element provider.

    Parameters
    ----------
    elType
        A string identifying the requested element formulation.
    provider
        The name of the element provider ot load.

    Returns
    -------
    type
        The element provider class type.
    """

    if provider is None:
        provider = "marmot"

    if strCaseCmp(provider, "edelweiss"):

        try:
            return eval(elLibrary[elType]["elClass"])
        except KeyError:
            raise Exception("Edelweiss element not found in library.")

    if provider.lower() == "marmot":
        try:
            from edelweissfe.elements.marmotelement.element import MarmotElementWrapper
        except ImportError:
            raise NotImplementedError("MarmotElementWrapper not found, maybe you built EdelweissFE without Marmot")

        return MarmotElementWrapper

    if provider.lower() == "marmotsingleqpelement":
        try:
            from edelweissfe.elements.marmotsingleqpelement.element import (
                MarmotMaterialWrappingElement,
            )
        except ImportError:
            raise NotImplementedError(
                "MarmotMaterialWrappingElement not found, maybe you built EdelweissFE without Marmot"
            )

        return MarmotMaterialWrappingElement

    else:
        raise Exception("This element provider doesn't exist!")
