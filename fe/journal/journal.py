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
Created on Sun Jan 15 11:30:59 2017

@author: Matthias Neuner
"""

import math


class Journal:
    """This class provides an interface to present messages to the user via console
    output and/or file output.
    Information messages can be sorted by the importance level.
    Suppressing certain levels of output is possible."""

    def __init__(self, verbose=True, outputFile=None, suppressFromLevel=3):
        self.suppressLvl = suppressFromLevel
        self.verbose = verbose
        self.setNewLineWidth(newWidth=100, leftColumn=80)

    def setNewLineWidth(self, newWidth=100, leftColumn=80):
        self.linewidth = newWidth
        self.leftColumn = leftColumn
        self.outputWidths = {}
        self.outputWidths[0] = leftColumn - 1
        self.outputWidths[1] = leftColumn - 4
        self.outputWidths[2] = leftColumn - 5
        self._rightColumn = self.linewidth - leftColumn

        self.errorMessageTemplate = " > > > {{:<{:}}}{{:>{:}}} < < < ".format(leftColumn - 12, self._rightColumn - 2)
        self.leveledOutput = {
            0: " {{:<{:}}}{{:>{:}}} ".format(leftColumn, self._rightColumn - 2),
            1: "   {{:<{:}}}{{:>{:}}} ".format(leftColumn - 2, self._rightColumn - 2),
            2: "     {{:<{:}}}{{:>{:}}} ".format(leftColumn - 4, self._rightColumn - 2),
        }

    def message(self, message, senderIdentification, level=1):
        while len(message) >= self.leftColumn:
            self.setNewLineWidth(self.linewidth + 5, self.leftColumn + 5)

        if level < self.suppressLvl:
            if self.verbose:
                print(self.leveledOutput[level].format(message, senderIdentification))

    def errorMessage(self, errorMessage, senderIdentification):
        print(self.errorMessageTemplate.format(errorMessage, senderIdentification))

    def printSeperationLine(
        self,
    ):
        if self.verbose:
            print("+" + "-" * (self.linewidth - 2) + "+")

    def printTable(self, table, senderIdentification, level=1, printHeaderRow=True):

        nCols = len(table[0])

        cellWidth = int(math.floor(self.outputWidths[level] / nCols - 1))
        rowBar = "+" + (("-" * cellWidth) + "+") * nCols
        rowString = ("|{:}".format("{:" + str(cellWidth) + "}")) * nCols + "|"

        if printHeaderRow:
            self.message(rowBar, senderIdentification, level)
        for row in table:
            self.message(rowString.format(*row), senderIdentification, level)

        self.message(rowBar, senderIdentification, level)

    def setVerbose(
        self,
    ):
        self.suppressLvl = 3

    def squelch(self, level):
        self.suppressLvl = level
