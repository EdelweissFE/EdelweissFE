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

from edelweissfe.models.femodel import FEModel
from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase

"""
Writes the (generated) mesh data to a file.

.. code-block:: console
    :caption: Example:

    *output, type=meshdatatofile, name=meshdata
        filename=myMesh.inc
"""

documentation = {
    "filename": "(optional), name of the file to write to; default: 'jobname'_mesh.inc",
}


class OutputManager(OutputManagerBase):
    """Simple status file writer for step, incrementation and iterations"""

    identification = "Meshdatatofile"
    printTemplate = "{:}, {:}: {:}"

    def __init__(self, name, model, fieldOutputController, journal, plotter):
        self.journal = journal
        self.filename = "{:}_mesh.inc".format(name)
        self.model = model

    def updateDefinition(self, **kwargs: dict):
        if "filename" in kwargs:
            self.filename = kwargs.get("filename")

        self.writeMeshDataToFile(self.model)

    def initializeJob(self):
        pass

    def initializeStep(self, step):
        pass

    def finalizeIncrement(self, statusInfoDict: dict = {}, **kwargs):
        pass

    def finalizeFailedIncrement(self, statusInfoDict: dict = {}):
        pass

    def finalizeStep(self):
        pass

    def finalizeJob(
        self,
    ):
        pass

    def writeMeshDataToFile(self, model: FEModel):
        """Write the mesh data to a file.

        Parameters
        ----------
        model
            The model dictionary containing the mesh data.
        """

        with open(self.filename, "w+") as f:
            # write nodes
            f.write("*NODE\n")
            for nodeID in model.nodes:
                f.write("{:},".format(nodeID))
                [f.write(" {:},".format(coord)) for coord in model.nodes[nodeID].coordinates]
                f.write("\n")

            # write elements
            # first, get all element types
            elementTypes = set()
            [elementTypes.add(element.elType) for element in model.elements.values()]

            for elementType in elementTypes:
                f.write("*ELEMENT, TYPE={:}\n".format(elementType))
                for elementID in model.elements:
                    if model.elements[elementID].elType == elementType:
                        f.write("{:>5},".format(elementID))
                        [f.write(" {:>5},".format(node.label)) for node in model.elements[elementID].nodes]
                        f.write("\n")

            # write node sets
            for nodeSetName in model.nodeSets:
                f.write("*NSET, NSET={:}\n".format(nodeSetName))
                counter = 0
                for node in model.nodeSets[nodeSetName]:
                    counter += 1
                    f.write(" {:>5},".format(node.label))
                    if counter % 16 == 0:
                        f.write("\n")
                f.write("\n")

            # write element sets
            for elementSetName in model.elementSets:
                f.write("*ELSET, ELSET={:}\n".format(elementSetName))
                counter = 0
                for element in model.elementSets[elementSetName]:
                    counter += 1
                    f.write(" {:>5},".format(element.elNumber))
                    if counter % 16 == 0:
                        f.write("\n")
                f.write("\n")

            # write surfaces
            for surfaceName in model.surfaces:
                f.write("*SURFACE, TYPE=ELEMENT, NAME={:}\n".format(surfaceName))
                for faceID in model.surfaces[surfaceName].keys():
                    elset = model.surfaces[surfaceName][faceID]
                    f.write("{elset}, {faceID}\n".format(elset=elset.name, faceID=faceID))
