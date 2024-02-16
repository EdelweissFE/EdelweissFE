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
# Created on Tue Aug 9 15:41:51 2022

"""
Find the node closest to a given spatial position, and store it in an existing or new node set.
"""

documentation = {"location": "The location of the node", "storeIn": "Store in this node set"}

from fe.points.node import Node
from fe.sets.nodeset import NodeSet
from fe.utils.misc import convertLinesToStringDictionary
from fe.journal.journal import Journal
from fe.models.femodel import FEModel

from fe.utils.exceptions import WrongDomain

import numpy as np

identification = "findclosestnode"


def generateModelData(generatorDefinition: dict, model: FEModel, journal: Journal):
    options = generatorDefinition["data"]
    options = convertLinesToStringDictionary(options)

    loc = np.fromstring(options["location"], sep=",", dtype=float)

    if len(loc) != model.domainSize:
        raise WrongDomain("Spatial dimension of specified location does not match model dimension")

    allNodes = np.asarray([n.coordinates for n in model.nodes.values()])

    differenceNorm = np.linalg.norm(allNodes - loc, axis=1)

    indexClosest = differenceNorm.argmin()

    closestNode = list(model.nodes.values())[indexClosest]

    if options["storeIn"] in model.nodeSets:
        raise Exception("Nodeset {options['storeIn']} already exists")

    model.nodeSets[options["storeIn"]] = NodeSet(options["storeIn"], [closestNode])

    return model
