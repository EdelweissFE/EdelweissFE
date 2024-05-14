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
#  Alexander Dummer alexander.dummer@uibk.ac.at
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

from collections import defaultdict
from typing import Union

from edelweissfe.config.timing import createTimingDict, timingTypes
from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase

"""
Prints the compute times per increment to the screen and writes them into a file (optional).

.. code-block:: console
    :caption: Example:

    *output, type=computetimemonitor, name=mycomputetimes
        export=myComputeTimes
"""

documentation = {
    "export": "(optional), filename if compute time should be written in a file",
}


class OutputManager(OutputManagerBase):
    identification = "ComputeTimeMonitor"

    def __init__(self, name, model, fieldOutputController, journal, plotter):
        self.journal = journal
        self.computingTimesOld = defaultdict(lambda: 0.0)
        self.stepcounter = 0
        self.exportFile = None

    def updateDefinition(self, **kwargs: dict):
        self.exportFile = kwargs.get("export")

        if self.exportFile is not None:
            with open(self.exportFile, "w+") as f:
                f.write("# \n# EdelweissFE: computing times per increment\n#\n")

    def initializeJob(self):
        pass

    def initializeStep(self, step):
        self.computingTimesOld = createTimingDict()
        self.stepcounter += 1

        if self.exportFile is not None:
            with open(self.exportFile, "a") as f:
                f.write("#\n# simulation step {:}\n#\n#".format(self.stepcounter))

                f.write(" {:<20} {:<20} {:<20}".format("increment", "simulation time", "inc compute time"))
                for timingType in timingTypes:
                    f.write(" {:<20}".format(timingType))
                f.write("\n#\n")

    def finalizeIncrement(
        self, statusInfoDict=defaultdict(lambda: 0.0), currentComputingTimes=createTimingDict(), **kwargs
    ):
        compTimeTotal, compTimeIndividual = self.computeIncrementComputingTimes(currentComputingTimes)

        self.printIncrementComputingTimes(compTimeTotal, compTimeIndividual)
        if self.exportFile is not None:
            self.writeIncrementComputingTimesToFile(compTimeTotal, compTimeIndividual, statusInfoDict)
        self.computingTimesOld = currentComputingTimes.copy()

    def finalizeFailedIncrement(
        self, statusInfoDict=defaultdict(lambda: 0.0), currentComputingTimes=defaultdict(lambda: 0.0), **kwargs
    ):
        compTimeTotal, compTimeIndividual = self.computeIncrementComputingTimes(currentComputingTimes)

        self.printIncrementComputingTimes(compTimeTotal, compTimeIndividual)
        if self.exportFile is not None:
            self.writeIncrementComputingTimesToFile(compTimeTotal, compTimeIndividual, statusInfoDict)
        self.computingTimesOld = currentComputingTimes.copy()

    def finalizeStep(
        self,
    ):
        pass

    def finalizeJob(
        self,
    ):
        pass

    def computeIncrementComputingTimes(self, currentComputingTimes: dict) -> Union[float, dict]:
        incrementComputingTimes = createTimingDict()
        incrementComputingTimeTotal = 0.0
        for key, val in currentComputingTimes.items():
            incrementComputingTimes[key] = val - self.computingTimesOld[key]
            incrementComputingTimeTotal += incrementComputingTimes[key]

        return incrementComputingTimeTotal, incrementComputingTimes

    def printIncrementComputingTimes(self, incCompTimeTotal, incCompTimesIndividual):
        self.journal.printTable(
            [
                (
                    "Time in {:}".format(k),
                    " {:.5e}s ({:>4.1f}%)".format(v, v / incCompTimeTotal * 100),
                )
                for k, v in incCompTimesIndividual.items()
            ],
            self.identification,
        )

    def writeIncrementComputingTimesToFile(self, incCompTimeTotal, incCompTimesIndividual, statusInfoDict):
        with open(self.exportFile, "a") as f:
            f.write("  {:<20}".format(int(statusInfoDict["inc"])))
            f.write(" {:<20.5e}".format(statusInfoDict["time end"]))
            f.write(" {:<20.5e}".format(incCompTimeTotal))
            for val in incCompTimesIndividual.values():
                f.write(" {:<20.5e}".format(val))
            f.write("\n")
