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
# Created on 2022-03-08

# @author: Paul Hofer
"""
Interface to Cubit. Generate mesh using Cubit .jou files.
"""

import os
import shlex

from edelweissfe.generators.abqmodelconstructor import AbqModelConstructor
from edelweissfe.utils.inputfileparser import parseInputFile
from edelweissfe.utils.misc import convertLinesToStringDictionary, strtobool

documentation = {
    "cubitCmd": "(Optional) Cubit executable; default=cubit",
    "jouFile": "Path to Cubit journal (.jou) file",
    "outFile": "(Optional) path to output file; default=mesh.inc",
    "vars": "(Optional) APREPRO variables as string '...' of comma-separated <key>=<value> pairs",
    "elType": "Specify element type for all sections",
    "elTypePerBlock": "Specify element type for each section as string '...' of comma-separated <key>=<value> pairs",
    "overwrite": "(Optional) overwrite existing outFiles; default=False",
    "runCubit": "(Optional) run Cubit GUI for debugging purposes; default=False",
    "silent": "(Optional) hide Cubit output; default=False",
}


def generateModelData(generatorDefinition, model, journal):
    options = generatorDefinition["data"]
    options = convertLinesToStringDictionary(options)

    # name = generatorDefinition.get("name", "cubit")

    cubitCmd = options.get("cubitCmd", "cubit")
    jouFile = options.get("jouFile")
    outFile = options.get("outFile", "mesh.inc")
    APREPROVars = options.get("APREPROVars")
    elType = options.get("elType")
    elTypePerBlock = options.get("elTypePerBlock")
    overwrite = options.get("overwrite", "True")
    runCubit = options.get("runCubit", "False")
    silent = options.get("silent", "False")

    # getElementClass(options["elType"], options.get("elProvider", None))

    generate = False
    if not os.path.exists(outFile) or strtobool(overwrite):
        generate = True

    if generate:
        cubitOptns = []
        cubitOptns.append("-information off")
        cubitOptns.append("-nojournal")

        if not strtobool(runCubit):
            cubitOptns.append("-batch")
            cubitOptns.append("-nographics")

        if APREPROVars:
            varStr = ""
            s = shlex.shlex(APREPROVars.replace(" ", ""), posix=True)
            s.whitespace_split = True
            s.whitespace = ","
            varDict = dict(item.split("=", 1) for item in s)

            for key, value in varDict.items():
                varStr += "{}={} ".format(key, value)
            cubitOptns.append(varStr)

        cubitOptns.append(jouFile)

        optnStr = " ".join(cubitOptns)
        cmd = " ".join([cubitCmd, optnStr])

        exportFile = "./exportAbaqus.jou"
        with open(exportFile, "w+") as f:
            f.write('export abaqus "{}" partial overwrite\n'.format(outFile))
        cmd = " ".join([cmd, exportFile])

        if strtobool(silent):
            cmd = " ".join([cmd, "> /dev/null"])

        os.system(cmd)
        os.remove(exportFile)

    fileDict = parseInputFile(outFile)

    if elType:
        for elDef in fileDict["*element"]:
            elDef["type"] = elType

    if elTypePerBlock:
        s = shlex.shlex(elTypePerBlock.replace(" ", ""), posix=True)
        s.whitespace_split = True
        s.whitespace = ","
        elDict = dict(item.split("=", 1) for item in s)
        for elDef in fileDict["*element"]:
            elSet = elDef["elset"]
            elDef["type"] = elDict[elSet]

    abqModelConstructor = AbqModelConstructor(journal)
    model = abqModelConstructor.createGeometryFromInputFile(model, fileDict)
    model = abqModelConstructor.createSectionsFromInputFile(model, fileDict)
    model = abqModelConstructor.createConstraintsFromInputFile(model, fileDict)

    return model
