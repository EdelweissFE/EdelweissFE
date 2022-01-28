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
"""
Created on Tue Jan  17 19:10:42 2017

@author: Matthias Neuner

Inputfileparser for inputfiles employing the Abaqus syntax.
"""
import numpy as np
from os.path import dirname, join
import textwrap
import shlex
from fe.utils.caseinsensitivedict import CaseInsensitiveDict

typeMappings = {
    "integer": int,
    "float": float,
    "string": lambda x: x,
    "numpy float array": lambda x: np.asarray(x, dtype=np.double),
    "numpy integer array": lambda x: np.asarray(x, dtype=np.int),
}

inputLanguage = {
    "*element": (
        "definition of element(s)",
        {
            "elSet": ("string", "name"),
            "type": ("string", "assign one of the types definied in the elementlibrary"),
            "provider": ("string", "provider (library) for the element type. Default: Marmot"),
            "data": ("numpy integer array", "Abaqus like element definiton lines"),
        },
    ),
    "*elSet": (
        "definition of an element set",
        {
            "elSet": ("string", "name"),
            "generate": ("string", "set True to generate from data line 1: start-element, end-element, step"),
            "data": ("string", "Abaqus like element set definiton lines"),
        },
    ),
    "*node": (
        "definition of nodes",
        {
            "nSet": ("string", "name"),
            "data": ("numpy float array", "Abaqus like node definiton lines: label, x, [y], [z]"),
        },
    ),
    "*nSet": (
        "definition of an element set",
        {
            "nSet": ("string", "name"),
            "generate": ("string", "set True to generate from data line 1: start-node, end-node, step"),
            "data": ("string", "Abaqus like node set definiton lines"),
        },
    ),
    "*surface": (
        "definition of surface set",
        {
            "type": ("string", "type of surface (currently 'element' only)"),
            "name": ("string", "name"),
            "data": ("string", "Abaqus like definition. Type 'element': elSet, faceID"),
        },
    ),
    "*section": (
        "definition of an section",
        {
            "name": ("string", "name"),
            "thickness": ("float", "associated element set"),
            "material": ("string", "associated id of defined material"),
            "data": ("string", "list of associated element sets"),
            "type": ("string", "type of the section"),
        },
    ),
    "*material": (
        "definition of a material",
        {
            "name": ("string", "name of the property"),
            "id": ("string", "name of the property"),
            "statevars": ("integer", "(deprecated and ignored) number of statevars"),
            "data": ("numpy float array", "material properties, multiline possible"),
        },
    ),
    "*fieldOutput": (
        "define fieldoutput, which is used by outputmanagers",
        {
            "jobName": ("string", "(optional), name of job, standard=defaultJob"),
            "data": ("string", "defintions lines for the output module"),
        },
    ),
    "*output": (
        "define an output module",
        {
            "name": ("string", "(optional), name of manager, standard=None"),
            "jobName": ("string", "(optional), name of job, standard=defaultJob"),
            "type": ("string", "output module "),
            "data": ("string", "defintions lines for the output module"),
        },
    ),
    "*job": (
        "definition of an analysis job",
        {
            "name": ("string", "(optional) name of job, standard = defaultJob"),
            "solver": ("string", "(optional) define solver, standard = NIST"),
            "linSolver": ("string", "(optional) define linear solver, standard = superlu"),
            "domain": ("string", "define spatial domain: 1d, 2d, 3d"),
            "startTime": ("float", "(optional) start time of job, standard = 0.0"),
        },
    ),
    "*step": (
        "definition of *job steps",
        {
            "stepLength": ("float", "time period of step"),
            "jobName": ("string", "(optional), name of job, standard=defaultJob"),
            "maxInc": ("float", "maximum size of increment"),
            "minInc": ("float", "minimum size of increment"),
            "maxNumInc": ("integer", "maximum number of increments"),
            "maxIter": ("integer", "maximum number of iterations"),
            "criticalIter": ("integer", "maximum number of iterations to prevent from increasing the increment"),
            "data": ("string", "define step actions, which are handled by the corresponding stepaction modules"),
        },
    ),
    "*updateConfiguration": (
        "update an configuration",
        {
            "configuration": ("string", " name of the modified settings category"),
            "data": ("string", "key=value pairs"),
        },
    ),
    "*modelGenerator": (
        "define a model generator, loaded from a module",
        {
            "generator": ("string", "generator module"),
            "name": ("string", "(optional) name of the generator"),
            "data": ("string", "key=value pairs"),
        },
    ),
    "*constraint": (
        "define a constraint",
        {
            "type": ("string", "constraint type"),
            "name": ("string", "(optional) name of the constraint"),
            "data": ("string", "definition of the constraint"),
        },
    ),
    "*configurePlots": (
        "customize the figures and axes",
        {
            "data": ("string", "key=value pairs for configuration of figures and axes"),
        },
    ),
    "*exportPlots": (
        "export your figures",
        {
            "data": ("string", "key=value pairs for exporting of figures and axes"),
        },
    ),
    "*include": (
        "(optional) load extra .inp file (fragment), use relative path to current .inp",
        {"input": ("string", "filename")},
    ),
    "*parameterIdentification": (
        "identify material parameter for given x and y data",
        {
            "xData": ("string", "filename where xData is given"),
            "yData": ("string", "filename where yData is given"),
            "plot": ("string", "True|False plot result with final parameters"),
            "id": ("string", "name of the property"),
            "jobName": ("string", "name of job"),
            "file": ("string", "filename where identified parameters are written"),
            "data": ("string", "key=value pairs"),
        },
    ),
}
inputLanguage_ = CaseInsensitiveDict()
for kw, (doc, opts) in inputLanguage.items():
    inputLanguage_[kw] = (doc, CaseInsensitiveDict(opts))
inputLanguage = inputLanguage_


def parseInputFile(fileName, currentKeyword=None, existingFileDict=None):
    """Parse an Abaqus like input file to generate an dictionary with its content"""

    if not existingFileDict:
        fileDict = CaseInsensitiveDict({key: [] for key in inputLanguage.keys()})
    else:
        fileDict = existingFileDict

    keyword = currentKeyword
    with open(fileName) as f:
        for l in f:
            lexer = shlex.shlex(l.strip(), posix=True)
            lexer.whitespace_split = True
            lexer.whitespace = ","
            lineElements = [x.strip() for x in lexer]

            if not lineElements or lineElements[0].startswith("**"):
                # line is a comment
                pass
            elif lineElements[0].startswith("*"):
                # line is keywordline
                lastkeyword = keyword
                keyword = lineElements[0]
                optionAssignments = lineElements[1:]

                objectentry = CaseInsensitiveDict()
                objectentry["data"] = []
                objectentry["inputFile"] = fileName  # save also the filename of the original inputfile!

                for ass in optionAssignments:
                    opts = ass.split("=")
                    optKey = opts[0].strip()
                    val = opts[1].strip()

                    doc, options = inputLanguage[keyword]
                    if optKey not in options:
                        raise KeyError('option "{:}" not valid for {:}'.format(optKey, keyword))

                    optionDataType, optionDoc = options[optKey]
                    try:
                        objectentry[optKey] = typeMappings[optionDataType](val)
                    except ValueError:
                        raise ValueError(
                            '{:}, option {:}: cannot convert "{:}" to {:}'.format(keyword, optKey, val, optionDataType)
                        )
                    except:
                        raise

                # special treatment for *include:
                if keyword == "*include":
                    includeFile = objectentry["input"]
                    parseInputFile(
                        join(dirname(fileName), includeFile), currentKeyword=lastkeyword, existingFileDict=fileDict
                    )
                    keyword = lastkeyword
                else:
                    fileDict[keyword].append(objectentry)

            else:
                # line is a dataline
                if "data" not in inputLanguage[keyword][1]:
                    raise KeyError("{:} expects no data lines".format(keyword))

                data = lineElements
                dataType = inputLanguage[keyword][1]["data"][0]

                try:
                    if "numpy" in dataType:
                        data = typeMappings[dataType](data)
                    else:
                        data = [typeMappings[dataType](d) for d in data]
                except ValueError:
                    raise ValueError("{:} data line: cannot convert {:} to {:}".format(keyword, data, dataType))

                fileDict[keyword][-1]["data"].append(data)
    return fileDict


def printKeywords():
    """print the input file language set"""
    kwString = "    {:}    "
    kwDataString = "        {:22}{:20}"

    wrapper = textwrap.TextWrapper(width=80, replace_whitespace=False)
    for kw, (kwDoc, optiondict) in sorted(inputLanguage.items()):
        wrapper.initial_indent = kwString.format(str(kw))
        wrapper.subsequent_indent = " " * len(wrapper.initial_indent)
        print(wrapper.fill(kwDoc))
        print("")
        for key in sorted(optiondict.keys()):
            optionName = key
            dType, description = optiondict[key]
            wrapper.initial_indent = kwDataString.format(str(optionName), dType)
            wrapper.subsequent_indent = " " * len(wrapper.initial_indent)
            print(wrapper.fill(description))
        print("\n")
