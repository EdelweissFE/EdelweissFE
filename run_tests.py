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

import os

import matplotlib

matplotlib.use("Agg")

from fe.fecore import finiteElementSimulation
from fe.utils.inputfileparser import parseInputFile
import numpy as np
import argparse
from timeit import default_timer as timer
from rich import print
import rich

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="validation script for FE analyses")
    parser.add_argument("--create", dest="create", action="store_true", help="create refernece solutions")
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

    testfilesDir = os.path.join(os.getcwd(), "testfiles")
    os.chdir(testfilesDir)
    testsDirs = next(os.walk("."))[1]

    if "all" not in tests:
        testsDirs = list(set(testsDirs).intersection(set(tests)))

    for directory in testsDirs:

        os.chdir(os.path.join(testfilesDir, directory))

        # no test.inp file is found
        if not os.path.exists(testfile):
            continue

        inputFile = parseInputFile(testfile)
        print("Test {:50}".format(directory), end="\r")

        try:
            tic = timer()
            success, U, P, fieldOutputController = finiteElementSimulation(inputFile, verbose=False, suppressPlots=True)
            toc = timer()

            if not args.create:
                if success:
                    UReference = np.loadtxt(referenceSolutionFile)
                    residual = U - UReference

                    if (np.max(np.abs(residual))) < 1e-6:
                        print("Test {:50} [green]PASSED[/] [{:2.1f}]".format(directory, toc - tic))
                    else:
                        print("Test {:50} [red]FAILED[/]".format(directory))
                    os.chdir("..")
                else:
                    print("Test {:50} [red]FAILED[/]: ".format(directory) + "Test not completed!")

            else:
                print("")
                np.savetxt(referenceSolutionFile, U)

        except Exception as e:
            print("Test {:50} [red]FAILED[/]: ".format(directory) + str(e))
            continue
