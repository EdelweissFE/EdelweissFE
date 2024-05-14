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
Created on Wed Jun 14 21:40:55 2017

@author: Matthias Neuner
"""

import argparse
import os
import sys
from timeit import default_timer as timer

import matplotlib
import numpy as np
from rich import print

from edelweissfe.drivers.inputfiledrivensimulation import finiteElementSimulation
from edelweissfe.utils.inputfileparser import parseInputFile

matplotlib.use("Agg")


def main():
    parser = argparse.ArgumentParser(description="validation script for FE analyses")
    parser.add_argument("testdirectory", help="The directory containing the testfiles")
    parser.add_argument(
        "--create",
        dest="create",
        action="store_true",
        help="create reference solutions",
    )
    parser.add_argument(
        "--tests",
        help="comma-separated list (without whitespace inbetween) with names of analyzed test files, "
        "e.g. MeshPlot,NodeForces, or simply type all. The names are case-sensitive.",
        type=str,
        default="all",
    )
    args = parser.parse_args()

    testfile = "test.inp"
    referenceSolutionFile = "U.ref"
    tests = [item for item in args.tests.split(",")]

    testfilesDir = os.path.abspath(args.testdirectory)
    os.chdir(testfilesDir)
    testsDirs = next(os.walk("."))[1]
    testsDirs = sorted(testsDirs, key=str.casefold)

    if "all" not in tests:
        testsDirs = list(set(testsDirs).intersection(set(tests)))

    failedTests = 0

    for directory in testsDirs:
        os.chdir(os.path.join(testfilesDir, directory))

        # no test.inp file is found
        if not os.path.exists(testfile):
            continue

        inputFile = parseInputFile(testfile)
        print("Test {:50}".format(directory), end="\r")

        try:
            tic = timer()
            model, fieldOutputController = finiteElementSimulation(inputFile, verbose=False, suppressPlots=True)
            toc = timer()

            u = [f["U"].flatten() for f in model.nodeFields.values()]
            sv = [v.value for v in model.scalarVariables.values()]

            U = np.hstack(u + sv).flatten()

            if not args.create:
                UReference = np.loadtxt(referenceSolutionFile)
                residual = U - UReference

                if (np.max(np.abs(residual))) < 1e-6:
                    print("Test {:50} [green]PASSED[/] [{:2.1f}]".format(directory, toc - tic))
                else:
                    print("Test {:50} [red]FAILED[/]".format(directory))
                    failedTests += 1
                os.chdir("..")
            else:
                print("")
                np.savetxt(referenceSolutionFile, U)

        except NotImplementedError as e:
            print("Test {:50} [grey]SKIPPED[/]: ".format(directory) + str(e))
            continue

        except Exception as e:
            print("Test {:50} [red]FAILED[/]: ".format(directory) + str(e))
            failedTests += 1
            continue

    print("[blue]Tests failed: {:3}[/]".format(failedTests))
    if failedTests > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
