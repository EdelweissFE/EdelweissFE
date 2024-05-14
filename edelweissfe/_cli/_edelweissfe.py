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
"""

import argparse
import sys
import warnings

import matplotlib
import numpy as np

from edelweissfe.drivers.inputfiledrivensimulation import finiteElementSimulation
from edelweissfe.utils.inputfileparser import (
    parseInputFile,
    printKeywords,
    printKeywordsRST,
)
from edelweissfe.utils.printdocumentation import printDocumentation

warnings.simplefilter("always", DeprecationWarning)


def main():
    parser = argparse.ArgumentParser(description="Batch computation for Edelweiss finite element jobs")

    parser.add_argument(
        "file",
        type=str,
        nargs="*",
    )
    parser.add_argument("--quiet", dest="verbose", action="store_false", help="suppress output")
    parser.add_argument("--noplot", dest="noplot", action="store_true", help="suppress plots")
    parser.add_argument(
        "--mplBackend",
        dest="mplBackend",
        default=None,
        type=str,
        help="define a matplotlib backend",
    )
    parser.add_argument(
        "--output",
        dest="output",
        default=None,
        type=str,
        help="write the final solution to a file",
    )
    parser.add_argument("--keywords", dest="kw", action="store_true", help="print keywords")
    parser.add_argument(
        "--keywordsRST",
        dest="kwRST",
        action="store_true",
        help="print keywords in RST format",
    )
    parser.add_argument("--doc=module", dest="doc", help="print keywords")
    args = parser.parse_args()

    fileList = args.file

    if not len(sys.argv) > 1:
        parser.print_help()
        exit(0)

    if args.kw:
        printKeywords()
        exit(0)

    if args.kwRST:
        printKeywordsRST()
        exit(0)

    if args.doc:
        printDocumentation(args.doc)
        exit(0)

    if args.mplBackend:
        print("Setting matplotlib backend to {:}".format(args.mplBackend))
        matplotlib.use(args.mplBackend)

    inputFiles = []

    for file in fileList:
        inputFiles.append(parseInputFile(file))

    for inputFile in inputFiles:
        model, fieldOutputController = finiteElementSimulation(
            inputFile, verbose=args.verbose, suppressPlots=args.noplot
        )

        if args.output:
            U = np.hstack(
                [f["U"].flatten() for f in model.nodeFields.values()]
                + [v.value for v in model.scalarVariables.values()]
            ).flatten()
            np.savetxt(args.output, U)


if __name__ == "__main__":
    main()
