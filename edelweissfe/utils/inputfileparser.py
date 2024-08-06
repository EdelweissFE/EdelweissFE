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
Inputfileparser for inputfiles employing an Abaqus-like syntax.
"""

import shlex
import textwrap
import warnings
from os.path import dirname, join

from edelweissfe.utils.caseinsensitivedict import CaseInsensitiveDict

typeMappings = {
    "integer": int,
    "float": float,
    "bool": bool,
    "string": lambda x: x,
    # "numpy float array": lambda x: np.asarray(x, dtype=np.double),
    # "numpy integer array": lambda x: np.asarray(x, dtype=np.int),
}


inputLanguage = {
    "*element": (
        "definition of element(s)",
        {
            "elSet": ("string", "name"),
            "type": ("string", "assign one of the types defined in the elementlibrary"),
            "provider": (
                "string",
                "provider (library) for the element type. Default: Marmot",
            ),
            "data": ("numpy integer array", "Abaqus like element definition lines"),
        },
    ),
    "*elSet": (
        "definition of an element set",
        {
            "elSet": ("string", "name"),
            "generate": (
                "string",
                "set True to generate from data line 1: start-element, end-element, step",
            ),
            "data": ("string", "Abaqus like element set definition lines"),
        },
    ),
    "*node": (
        "definition of nodes",
        {
            "nSet": ("string", "name"),
            "data": (
                "numpy float array",
                "Abaqus like node definition lines: label, x, [y], [z]",
            ),
        },
    ),
    "*nSet": (
        "definition of an element set",
        {
            "nSet": ("string", "name"),
            "generate": (
                "string",
                "set True to generate from data line 1: start-node, end-node, step",
            ),
            "data": ("string", "Abaqus like node set definition lines"),
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
        "definition of a section",
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
            "provider": ("string", "material provider"),
        },
    ),
    "*fieldOutput": (
        "define fieldoutput, which is used by outputmanagers",
        {
            "data": ("string", "definition lines for the output module"),
        },
    ),
    "*analyticalField": (
        "define an analytical field",
        {
            "name": ("string", "name of analytical field"),
            "type": (
                "string",
                "type of analytical field (currently 'expression' only)",
            ),
            "data": ("string", "definition"),
        },
    ),
    "*output": (
        "define an output module",
        {
            "name": ("string", "(optional), name of manager, standard=None"),
            "type": ("string", "output module "),
            "data": ("string", "definition lines for the output module"),
        },
    ),
    "*job": (
        "definition of an analysis job",
        {
            "name": ("string", "(optional) name of job, standard = defaultJob"),
            "domain": ("string", "define spatial domain: 1d, 2d, 3d"),
            "solver": ("string", "(deprecated) define the solver to be used"),
            "startTime": ("float", "(optional) start time of job, standard = 0.0"),
        },
    ),
    "*solver": (
        "definition of solver",
        {
            "name": ("string", "Name of this solver"),
            "solver": ("string", "Solvertype"),
            "data": (
                "string",
                "define options which are passed to the respective solver instance.",
            ),
        },
    ),
    "*step": (
        "definition of job steps",
        {
            "stepLength": ("float", "time period of step"),
            "startInc": ("float", "size of the start increment"),
            "maxInc": ("float", "maximum size of increment"),
            "minInc": ("float", "minimum size of increment"),
            "maxNumInc": ("integer", "maximum number of increments"),
            "maxIter": ("integer", "maximum number of iterations"),
            "type": ("string", "(optional) define step type, default = AdaptiveStep"),
            "solver": ("string", "(optional) solver to be used"),
            "criticalIter": (
                "integer",
                "maximum number of iterations to prevent from increasing the increment",
            ),
            "data": (
                "string",
                "define step actions, which are handled by the corresponding stepaction modules",
            ),
        },
    ),
    "*updateConfiguration": (
        "update a configuration",
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
            "executeAfterManualGeneration": (
                "bool",
                "(optional) Delay the execution of the generator after the manual creation of the mesh, default=False",
            ),
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
}

inputLanguage_ = CaseInsensitiveDict()
for kw, (doc, opts) in inputLanguage.items():
    inputLanguage_[kw] = (doc, CaseInsensitiveDict(opts))
inputLanguage = inputLanguage_


def parseInputFile(
    fileName: str,
    currentKeyword: str = None,
    existingFileDict: CaseInsensitiveDict = None,
) -> CaseInsensitiveDict:
    """Parse an Abaqus-like input file to generate a dictionary with its content.

    Parameters
    ----------
    fileName
        The name of the file to parse.
    currentKeyword
        If nested parsing is performed by using ``*include``, this option tells which
        keyword is currently active.
    existingFileDict
        An existing dictionary to append. If Nonde, a new dictionary is created.

    Returns
    -------
    CaseInsensitiveDict
        The parsed input file.
    """

    if not existingFileDict:
        fileDict = CaseInsensitiveDict({key: [] for key in inputLanguage.keys()})
    else:
        fileDict = existingFileDict

    keyword = currentKeyword
    with open(fileName) as f:
        # filter out empty lines and comments
        lines = (line.strip() for line in f)
        lines = (line for line in lines if line and not line.startswith("**"))

        for line in lines:
            if line.startswith("*"):
                lexer = shlex.shlex(line, posix=True)
                lexer.whitespace_split = True
                lexer.whitespace = ","

                lineElements = [x.strip() for x in lexer]

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
                    except Exception as e:
                        raise e

                # special treatment for *include:
                if keyword.lower() == "*include":
                    includeFile = objectentry["input"]
                    parseInputFile(
                        join(dirname(fileName), includeFile),
                        currentKeyword=lastkeyword,
                        existingFileDict=fileDict,
                    )
                    keyword = lastkeyword
                else:
                    fileDict[keyword].append(objectentry)

            else:
                # line is a dataline
                if "data" not in inputLanguage[keyword][1]:
                    raise KeyError("{:} expects no data lines".format(keyword))

                fileDict[keyword][-1]["data"].append(line)

    # raise deprecation warning if deprecated jobName option is set in keywords
    keywords = ["*step", "*fieldOutput", "*output"]
    for keyword in keywords:
        for entry in fileDict[keyword]:
            if "jobName" in entry:
                warnings.warn(
                    f'Option "jobName" in {keyword} is deprecated and will be removed in future',
                    DeprecationWarning,
                    stacklevel=2,
                )

    for entry in fileDict["*job"]:
        if "solver" in entry:
            warnings.warn(
                "Warning, defining a Solver in *job is deprecated; Define solver using *solver keyword",
                DeprecationWarning,
                stacklevel=2,
            )

    return fileDict


def printKeywords():
    """Print the input file language set."""

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


def printKeywordsRST():
    """Print the input file language set in an RST conform format."""

    for kw, (kwDoc, optiondict) in sorted(inputLanguage.items()):
        print(".. list-table:: " + "``{:}`` : {:}".format(kw, kwDoc))
        print("    :width: 100%")
        print("    :widths: 25 25 40")
        print("    :header-rows: 1")
        print(" ")
        print("    * - Option")
        print("      - Type")
        print("      - Description")
        for key in sorted(optiondict.keys()):
            optionName = key
            dType, description = optiondict[key]

            print("    * - ``{:}``".format(optionName))
            print("      - ``{:}``".format(dType))
            print("      - " + description)
        print(" ")
