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

from time import time
from collections import defaultdict


class _PerformanceTimerBranch(defaultdict):
    def __init__(self):
        self.time = float()  #: the measured time for this branch.
        self._tic = None

        super().__init__(_PerformanceTimerBranch)

    def tic(
        self,
    ):
        """
        Start measuring time.
        """
        self._tic = time()

    def toc(
        self,
    ):
        """
        Stop measuring time.
        """
        self.time += time() - self._tic


times = _PerformanceTimerBranch()
"""The global dictionary of measured computations times."""


def timeit(category: str, *subCategories: list[str]):
    """Decorator function for indicating that a function should be timed.

    Parameters
    ----------
    category
        The top level category for storing the measured time.
    *args
        The optional list of sub categories.
    """

    def outer(theFunction):
        from time import time

        timer = times[category]
        for sub in subCategories:
            timer = timer[sub]

        def inner(*args, **kwargs):
            timer.tic()
            try:
                return theFunction(*args, **kwargs)
            finally:
                timer.toc()

        return inner

    return outer
