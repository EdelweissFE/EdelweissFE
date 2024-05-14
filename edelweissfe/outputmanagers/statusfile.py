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

from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase

"""
Writes a status file during the analysis.

.. code-block:: console
    :caption: Example:

    *output, type=statusfile, name=status
        filename=myStatus.sta
"""

documentation = {
    "filename": "(optional), custom filename for the status file; default: 'jobname'.sta",
}


class OutputManager(OutputManagerBase):
    """Simple status file writer for step, incrementation and iterations"""

    identification = "Statusfile"
    printTemplate = "{:}, {:}: {:}"

    def __init__(self, name, model, fieldOutputController, journal, plotter):
        self.journal = journal
        # self.filename = "{:}.sta".format(jobInfo.get("name", jobInfo["inputfile"].rstrip(".inp")))
        self.filename = "job.sta"
        self.statusFileExists = False

    def updateDefinition(self, **kwargs: dict):
        if "filename" in kwargs.keys():
            self.filename = kwargs.get("filename")

        self.statusFileExists = False

    # def initializeSimulation(self, model):
    #     pass

    def initializeJob(self):
        pass

    def initializeStep(self, step):
        pass

    def finalizeIncrement(self, statusInfoDict: dict = {}, **kwargs):
        self.writeStatusFile(statusInfoDict)

    def finalizeFailedIncrement(self, statusInfoDict: dict = {}, **kwargs):
        self.writeStatusFile(statusInfoDict)

    def finalizeStep(self):
        pass

    def finalizeJob(self):
        pass

    def writeStatusFile(self, statusInfoDict):
        """Write the status to a file.

        Parameters
        ----------
        statusInfoDict
            A dictionary containing information about the simulation status.
        """
        d = statusInfoDict

        if not self.statusFileExists:
            with open(self.filename, "w+") as f:
                f.write("#\n")
                f.write("# This is a status file of EdelweissFE.\n")
                f.write("#\n")
                f.write(
                    "#{: >5}{: >6}{: >6}{: >10}{: >12}{: >12}    {:<}\n".format(
                        "step",
                        "inc",
                        "iters",
                        "converged",
                        "time inc",
                        "time end",
                        "notes",
                    )
                )
                f.write("#\n")
            self.statusFileExists = True

        with open(self.filename, "a") as f:
            f.write(
                "{: >6}{: >6}{: >6}{: >10}{: >12.3e}{: >12.3e}    # {:<s}\n".format(
                    d["step"],
                    d["inc"],
                    d["iters"],
                    d["converged"],
                    d["time inc"],
                    d["time end"],
                    d["notes"],
                )
            )
