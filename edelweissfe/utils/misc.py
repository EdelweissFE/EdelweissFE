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
Created on Mon Apr 18 17:36:07 2016

@author: Matthias Neuner
"""

import shlex

import numpy as np


def flagDict(configLine):
    parts = [x.strip() for x in configLine.split("=")]
    opt = parts[0]
    val = True if (len(parts) > 1 and parts[1] == "True") else False
    return {opt: val}


def splitLineAtCommas(line: str) -> list:
    """Split a line at commas and strip the individual parts.

    Parameters
    ----------
    line
        The line to be split.

    Returns
    -------
    list
        The list of parts.
    """

    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.whitespace = ","

    lineElements = [x.strip() for x in lexer]

    return lineElements


def convertAssignmentsToStringDictionary(assignments: list) -> dict:
    """Create a dictionary from a list of assignments in
    the form a=b.

    Parameters
    ----------
    assignments
        The list of assignments.

    Returns
    -------
    dict
        The resulting dictionary.
    """

    # resultDict = CaseInsensitiveDict()
    resultDict = dict()
    for entry in assignments:
        parts = [x.strip() for x in entry.split("=")]
        opt = parts[0]
        val = "=".join(parts[1:]) if len(parts) > 1 else "True"
        resultDict[opt] = val

    return resultDict


def convertLinesToMixedDictionary(lines: list) -> dict:
    """Create a mixed dictionary from a list of strings containing multiple assignments in
    the form a=b. All strings containing evaluatable values will be evaluated.
    All other dictionary values stay the same.

    Parameters
    ----------
    lines
        The list of strings containing the assignments.

    Returns
    -------
    dict
        The resulting dictionary.
    """

    dictionary = convertLineToStringDictionary(",".join(lines)).copy()
    for key, value in dictionary.items():
        try:
            dictionary[key] = eval(value)
        except NameError:
            pass

    return dictionary


def convertLineToStringDictionary(line: str) -> dict:
    """Create a dictionary from a string containing multiple assignments in
    the form a=b.

    Parameters
    ----------
    line
        The string containing the assignments.

    Returns
    -------
    dict
        The resulting dictionary.
    """

    lineElements = splitLineAtCommas(line)

    return convertAssignmentsToStringDictionary(lineElements)


def convertLinesToStringDictionary(lines: list) -> dict:
    """Create a dictionary from a list of strings containing multiple assignments in
    the form a=b.

    Parameters
    ----------
    lines
        The list of strings containing the assignments.

    Returns
    -------
    dict
        The resulting dictionary.
    """

    return convertLineToStringDictionary(",".join(lines))


def convertLinesToFlatArray(lines: list, dtype: type = float) -> np.ndarray:
    """Create a 1D numpy array from a list of lines with elements separated by commas.

    Parameters
    ----------
    lines
        The list of strings containing the elements.

    Returns
    -------
    np.ndarray
        The resulting 1D array.
    """

    theLines = [np.asarray(splitLineAtCommas(line), dtype=dtype) for line in lines]
    return np.hstack(theLines)


def strCaseCmp(str1, str2):
    return str1.casefold() == str2.casefold()


def strToSlice(string):
    if ":" in string:
        idcs = [int(i) for i in string.split(":")]
        a, b = idcs
        return slice(a, b)
    else:
        return slice(int(string), int(string) + 1)


def strToRange(string):
    if ":" in string:
        idcs = [int(i) for i in string.split(":")]
        a, b = idcs
        return range(a, b)
    else:
        return range(int(string))


def isInteger(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def filterByJobName(canditates, jobName):
    return [cand for cand in canditates if "jobName" not in cand or cand["jobName"] == jobName]


def mergeNumpyDataLines(multiLineData: np.ndarray) -> np.ndarray:
    """Flatten a numpy array."""
    flattenedMatProps = [p for row in multiLineData for p in row]
    return np.array(flattenedMatProps, dtype=float)


def strtobool(val: str) -> bool:
    """-- Implementation from deprecated module distutils.utils --
    Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.

    Parameters
    ----------
    val
        The string representing the truth value:

    Returns
    -------
    bool
        The truth value.
    """
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return 1
    elif val in ("n", "no", "f", "false", "off", "0"):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))
