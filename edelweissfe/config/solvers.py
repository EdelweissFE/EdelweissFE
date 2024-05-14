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
# Created on Fri Feb 10 19:20:25 2017

# @author: Matthias Neuner
"""
Currently, EdelweissFE provides

 * a nonlinear implicit static solver (NIST),
 * a parallel nonlinear implicit static solver (NISTParallel),
 * a parallel nonlinear implicit static solver tuned for marmot elements (NISTParallelForMarmotElements),
 * and a parallel arc length solver (NISTPArcLength).

Choose the solver in the ``*solver`` definition:

.. code-block:: edelweiss

    *solver, name=mySolver, solver=NISTParallel
"""

import importlib

solverLibrary = {
    "NIST": "nonlinearimplicitstatic",
    "NISTParallel": "nonlinearimplicitstaticparallelmk2",
    "NISTParallelForMarmotElements": "nonlinearimplicitstaticparallel",
    "NISTPArcLength": "nonlinearimplicitstaticparallelarclength",
}


def getSolverByName(name: str) -> type:
    """Get the class type of the requested solver.

    Parameters
    ----------
    name
        The name of the solver to load.

    Returns
    -------
    type
        The solver class type.
    """

    try:
        solverType = solverLibrary[name]
    except KeyError:
        raise KeyError(f"Solver {name} not found in library. Available solvers: " + ", ".join(solverLibrary.keys()))

    solver = importlib.import_module("edelweissfe.solvers.{:}".format(solverType))
    return getattr(solver, name)
