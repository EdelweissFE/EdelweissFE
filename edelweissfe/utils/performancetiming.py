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

import inspect
from collections import defaultdict
from time import time

from prettytable import PrettyTable


class _PerformanceTimerBranch(defaultdict):
    def __init__(self):
        self.time = float()  #: the measured time for this branch.
        self.calls = int()
        self._tic = None

        super().__init__(_PerformanceTimerBranch)

    def tic(
        self,
    ):
        """
        Start measuring time.
        """
        self._tic = time()
        self.calls += 1

    def toc(
        self,
    ):
        """
        Stop measuring time.
        """
        self.time += time() - self._tic


times = _PerformanceTimerBranch()
"""The global dictionary of measured computations times."""


class timeit:
    """Decorator class for performance timing of functions.
    This decorator has a runtime memory, i.e., it is aware of the stack level
    of nested timed functions.

    Parameters
    ----------
    category
        The category for storing the measured time.
    """

    _currentStackLevel = times

    def __init__(self, category: str):
        self._category = category
        self._parentStackLevel = None

    def __call__(self, theFunction):
        def wrapper(*args, **kwargs):
            self._parentStackLevel = timeit._currentStackLevel
            timer = timeit._currentStackLevel[self._category]
            timeit._currentStackLevel = timer

            timer.tic()
            try:
                return theFunction(*args, **kwargs)
            finally:
                timer.toc()
                timeit._currentStackLevel = self._parentStackLevel

        wrapper.__doc__ = theFunction.__doc__
        wrapper.__module__ = theFunction.__module__
        wrapper.__signature__ = inspect.signature(theFunction)

        return wrapper


def _makeTable(branch: _PerformanceTimerBranch, level: int, maxLevels: int) -> list[tuple]:
    """Recursive function for creating a table of the measured times.

    Parameters
    ----------
    branch
        The current active branch.
    levels
        The current level.
    maxLevels
        The maximum number of stack levels considered in the table.

    Returns
    -------
    list[tuple]
        The table in list format containing columns as tuples."""

    table = []
    for k, v in branch.items():
        table.append((level, k, v.time, v.calls))
        if level < maxLevels and len(v):
            table += _makeTable(v, level + 1, maxLevels)

    return table


def makePrettyTable(maxLevels: int = 4) -> PrettyTable:
    """Create a pretty formatted table of the measured times.

    Parameters
    ----------
    maxLevels
        The maximum number of stack levels considered in the table.

    Returns
    -------
    PrettyTable
        The table in pretty format."""

    theTable = _makeTable(times, 0, maxLevels)

    prettytable = PrettyTable()
    prettytable.field_names = ["function", "acc. runtime", "calls"]
    prettytable.align = "l"

    for level, cat, t, calls in theTable:
        prettytable.add_row(
            (
                "{:}{:}".format(" " * level, cat),
                "{:}{:10.4f}s".format(" " * level, t),
                calls,
            )
        )

    return prettytable
