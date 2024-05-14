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
Created on Sat Feb 10 10:27:25 2018

@author: Matthias Neuner
"""


def getDefaultLinSolver():
    try:
        from edelweissfe.linsolve.pardiso.pardiso import pardisoSolve

        return pardisoSolve
    except ImportError:
        from scipy.sparse.linalg import spsolve

        return lambda A, b: spsolve(A, b, use_umfpack=False)


def getLinSolverByName(linsolverName):
    if linsolverName.lower() == "superlu":
        from scipy.sparse.linalg import spsolve

        return lambda A, b: spsolve(A, b, use_umfpack=False)
    elif linsolverName.lower() == "umfpack":
        from scipy.sparse.linalg import spsolve

        return lambda A, b: spsolve(A, b, use_umfpack=True)
    elif linsolverName.lower() == "pardiso":
        from edelweissfe.linsolve.pardiso.pardiso import pardisoSolve

        return pardisoSolve
    elif linsolverName.lower() == "panuapardiso":
        from edelweissfe.linsolve.panuapardiso.panuapardiso import panuaPardisoSolve

        return panuaPardisoSolve
    elif linsolverName.lower() == "amgcl":
        from edelweissfe.linsolve.amgcl.amgcl import amgclSolve

        return amgclSolve
    elif linsolverName.lower() == "klu":
        from edelweissfe.linsolve.klu.klu import kluSolve

        return kluSolve
    elif linsolverName.lower() == "petsclu":
        from edelweissfe.linsolve.petsclu.petsclu import petscluSolve

        return petscluSolve
    elif linsolverName.lower() == "mumps":
        from edelweissfe.linsolve.mumps.mumps import mumpsSolve

        return mumpsSolve
    else:
        raise AttributeError("invalid linear solver {:} requested".format(linsolverName))
