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
# Created on Sun Jan 15 11:30:59 2017

# @author: Matthias Neuner
"""
This module provides journaling capabilities,
which can be used troughout EdelweissFE.
"""

import math

from prettytable import PrettyTable


class Journal:
    """This class provides an interface to present messages to the user via console
    output and/or file output.
    Information messages can be sorted by the importance level.
    Suppressing certain levels of output is possible.

    Parameters
    ----------
    verbose
        Print to the terminal.
    outputFile
        Write to a file.
    suppressFromLevel
        Suppress certain levels.
    """

    def __init__(
        self,
        verbose: bool = True,
        outputFile: str = None,
        suppressFromLevel: int = 3,
    ):
        self._suppressLvl = suppressFromLevel
        self._verbose = verbose
        self._fileHandle = None
        self.setNewLineWidth(newWidth=100, leftColumn=80)

    def setFileOutput(self, fileHandle=None):
        if fileHandle:
            self._fileHandle = fileHandle
        else:
            self._fileHandle = None

    def setNewLineWidth(self, newWidth: int = 100, leftColumn: int = 80):
        """Set the line width of the log file.

        Parameters
        ----------
        newWidth
            New total width.
        leftColumn
            Width of the main column.
        """

        self.linewidth = newWidth
        self.leftColumn = leftColumn
        self.leftColumnMaxSize = 120
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

    def message(self, message: str, senderIdentification: str, level: int = 1):
        """Write message to log.

        Parameters
        ----------
        message
            The message.
        senderIdentification
            The name of the sender.
        level
            Level of message.
        """

        lines = message.splitlines()
        for line in lines:
            while len(line) > self.leftColumn and self.leftColumn < self.leftColumnMaxSize:
                self.setNewLineWidth(self.linewidth + 5, self.leftColumn + 5)

            theLine = self.leveledOutput[level].format(line, senderIdentification)
            if level < self._suppressLvl:
                if self._verbose:
                    print(theLine)

            if self._fileHandle:
                self._fileHandle.write(theLine + "\n")

    def errorMessage(self, errorMessage: str, senderIdentification: str):
        """Print an error message.

        Parameters
        ----------
        errorMessage
            The message.
        senderIdentification
            The name of the sender.
        """

        theLine = self.errorMessageTemplate.format(errorMessage, senderIdentification)
        print(theLine)

        if self._fileHandle:
            self._fileHandle.write(theLine + "\n")

    def printSeperationLine(
        self,
    ):
        """Write a seperation file to log."""

        theLine = "+" + "-" * (self.linewidth - 2) + "+"
        if self._verbose:
            print(theLine)

        if self._fileHandle:
            self._fileHandle.write(theLine + "\n")

    def printTable(
        self,
        table: list[list[str]],
        senderIdentification: str,
        level: int = 1,
        printHeaderRow: bool = True,
    ):
        """Print a pretty table.

        Parameters
        ----------
        table
            The table as nested list.
        senderIdentification
            The name of the sender.
        level
            The level of the mesage.
        printHeaderRow
            Special formatting for first row.
        """

        nCols = len(table[0])

        cellWidth = int(math.floor(self.outputWidths[level] / nCols)) - 3
        rowBar = "+" + (("-" * cellWidth) + "+") * nCols
        rowString = ("|{:}".format("{:" + str(cellWidth) + "}")) * nCols + "|"

        if printHeaderRow:
            self.message(rowBar, senderIdentification, level)
        for row in table:
            self.message(rowString.format(*row), senderIdentification, level)

        self.message(rowBar, senderIdentification, level)

    def printPrettyTable(self, prettyTable: PrettyTable, senderIdentification: str):
        prettyTable.min_table_width = self.leftColumn
        prettyTable.max_table_width = self.leftColumn
        self.message(str(prettyTable), senderIdentification, level=0)

    def setVerbose(
        self,
    ):
        """Set highest verbosity."""

        self._suppressLvl = 3

    def squelch(self, level):
        """Suppress all messages.

        Parameters
        ----------
        level
            The priority level of the message.
        """

        self._suppressLvl = level
