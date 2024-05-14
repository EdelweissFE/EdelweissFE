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
# Created on Mon Jul 24 15:30:36 2017

# @author: Matthias Neuner

import itertools
import os
import sys
from distutils.spawn import find_executable

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from edelweissfe.journal.journal import Journal

"""
Plotting in EdelweissFE through matplotlib can be performed easily
through a global plotting instance, which is passed to all output managers.
Output managers may use the plotter for visualizing all kinds of outputs.

The style of the plots can be configured through
the keyword ``*configurePlots``:

.. code-block:: edelweiss

    *configurePlots,
        figure=1, axSpec=211,
        figure=3, axSpec=111,  xLabel=U1, yLabel=P2, flipX=True
        figure=3, axSpec=111,  xLabel=U1, yLabel=P2, flipX=True
        figure=4, axSpec=111, aspect=equal

Plots can be exported to .pdf (and png) files using the ``*exportPlots`` keyword:

.. code-block:: edelweiss

    *exportPlots,
        figure=1, fileName=fig1, png=True
        figure=2, fileName=fig2

The default style of plots can be configured via a classical rcparams.py file,
which may be located in the current working directory.
"""

documentation_configurePlots = {
    "figure": "The figure to be configured",
    "xLimits": "Specify x axis limits",
    "yLimits": "Specify y axis limits",
    "xLabel": "Specify x axis label",
    "yLabel": "Specify y axis label",
    "flipX": "Flip x axis",
    "flipY": "Flip y axis",
    "aspect": "Set aspect",
    "grid": "Switch grid",
}

documentation_exportPlots = {
    "figure": "The figure to be exported",
    "fileName": "The export file name",
    "width": "Width of the figure",
    "heightRatio": "Ratio of height/width",
    "png": "Set true to export .png additionally",
    "scale": "Scale the figure",
}

defaultMarkerCycle = itertools.cycle(("o", "v", "D", "s", "^"))
defaultLinePlotColorCycle = itertools.cycle(("k"))
defaultScatterPlotColorCycle = itertools.cycle(("k", "r", "b", "g"))
defaultLineStyleCycle = itertools.cycle(("-", "0-3-2", "0-3-1-1-1", "0-1-1"))


class Plotter:
    """
    The unified plotter, which can be accessed and used by all outputmanagers.

    Parameters
    ----------
    journal
        The journal instance for logging.
    plotConfigurations
        The list of dictionaries configuring individual plots.
    exportJobs
        The list of jobs to export plots at the end of a simulations.
    """

    def __init__(self, journal: Journal, plotConfigurations: list[dict], exportJobs: list[dict]):
        self.plotConfigurations = plotConfigurations
        self.exportJobs = exportJobs

        latexAvailable = False
        if find_executable("latex"):
            latexAvailable = True

        self.rcParams = {
            "pgf.texsystem": "pdflatex",  # change this if using xetex or lautex
            "text.usetex": latexAvailable,  # use LaTeX to write all text
            "text.latex.preamble": r" \usepackage[utf8]{inputenc} \usepackage{amsmath} \usepackage{amssymb} \usepackage{mathpazo} \usepackage{siunitx}",  # noqa: E501
            "font.family": "serif",
            #                "font.serif": [],  # blank entries should cause plots to inherit fonts from the document
            #                "font.sans-serif": [],
            "font.monospace": [],
            "axes.labelsize": 10,  # LaTeX default is 10pt font.
            "font.size": 10,
            "legend.fontsize": 8,  # Make the legend/label fonts a little smaller
            "legend.numpoints": 1,
            "legend.labelspacing": 0.2,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "lines.linewidth": 1,
            "lines.markeredgewidth": 0.4,
            "lines.markersize": 4,
            "axes.unicode_minus": False,
        }
        plt.close("all")

        if os.path.isfile("plotterRcParams.py"):
            print("found plotterRcParams.py")
            sys.path.append(os.getcwd())
            import plotterRcParams

            self.rcParams.update(plotterRcParams.rcParams)

        self.figsWithAxes = {}  # {figID : (figure, {axesDict})}
        matplotlib.rcParams.update(self.rcParams)

    def getAx(self, figureID: int = 0, axSpec: int = 111) -> matplotlib.axes.Axes:
        """Get or create a figure with axes if it doesn't exist yet.

        Parameters
        ----------
        figureID
            The matplotlib figure ID.
        axSpec
            The matplotlib axis specification.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib Axes instance.
        """

        if figureID not in self.figsWithAxes:
            self.figsWithAxes[figureID] = (plt.figure(figureID), {})

        fig, axes = self.figsWithAxes[figureID]
        axSpec = axSpec
        if axSpec not in axes:
            self.figsWithAxes[figureID][1][axSpec] = fig.add_subplot(int(axSpec))
            self.figsWithAxes[figureID][1][axSpec].grid(True)

        return self.figsWithAxes[figureID][1][axSpec]

    def getFig(self, figureID: int = 0) -> matplotlib.figure.Figure:
        """Get or create a figure doesn't exist yet.

        Parameters
        ----------
        figureID
            The matplotlib figure ID.
        Returns
        -------
        matplotlib.figure.Figure
            The matplotlib Figure instance.
        """

        if figureID not in self.figsWithAxes:
            self.figsWithAxes[figureID] = (plt.figure(figureID), {})

        return self.figsWithAxes[figureID][0]

    def plotXYData(
        self,
        x: np.ndarray,
        y: np.ndarray,
        figureID: int = 1,
        axSpec: int = 111,
        plotOptions: dict = None,
    ):
        """Plots a single curve.

        Parameters
        ----------
        x
            The x data.
        y
            The y data.
        figureID
            The matplotlib figure ID.
        axSpec
            The matpotlib axes.
        plotOptions
            A dictionary with additional plot options in matplotlib format.
        """

        ax = self.getAx(figureID, axSpec)

        # flatten arrays
        x = np.asarray(x).flatten()
        y = np.asarray(y).flatten()

        plotDefinition = {}
        if "label" in plotOptions:
            plotDefinition["label"] = plotOptions["label"]
        if plotOptions.get("plotType", "linePlot") == "scatter":
            plotDefinition["marker"] = plotOptions.get("marker", None) or next(defaultMarkerCycle)
            if "ms" in plotOptions:
                plotDefinition["markersize"] = plotOptions.get("ms", None)
            plotDefinition["markerfacecolor"] = plotOptions.get("markerfacecolor", "None")
            plotDefinition["markeredgecolor"] = plotOptions.get("c", False) or next(defaultScatterPlotColorCycle)
            plotDefinition["linestyle"] = ""
        else:
            plotDefinition["c"] = plotOptions.get("c", False) or next(defaultLinePlotColorCycle)
            plotDefinition["linestyle"] = plotOptions.get("ls", False) or next(defaultLineStyleCycle)
            lsParts = plotDefinition["linestyle"].split("-")
            if len(lsParts) >= 3:
                plotDefinition["linestyle"] = (
                    float(lsParts[0]),
                    [int(onOff) for onOff in lsParts[1:]],
                )
        ax.plot(x, y, **plotDefinition)

    def configurePlotter(self):
        """Set global options of the plotter."""

        for configEntry in self.plotConfigurations:
            ax = self.figsWithAxes[(configEntry["figure"])][1][configEntry["axSpec"]]

            if "xLimits" in configEntry:
                limits = [float(x) for x in configEntry["xLimits"].split("_")]
                ax.set_xlim(limits)

            if "yLimits" in configEntry:
                limits = [float(x) for x in configEntry["yLimits"].split("_")]
                ax.set_ylim(limits)

            if "xLabel" in configEntry:
                ax.set_xlabel(configEntry["xLabel"])
            if "yLabel" in configEntry:
                ax.set_ylabel(configEntry["yLabel"])
            if "flipX" in configEntry:
                ax.invert_xaxis()
            if "flipY" in configEntry:
                ax.invert_yaxis()
            if "aspect" in configEntry:
                ax.set_aspect(configEntry["aspect"])
            if "grid" in configEntry:
                ax.grid()

    def exportPlots(self):
        """Export all plots according to the export job definitions."""

        for exportJob in self.exportJobs:
            self.exportFigure(
                exportJob.get("fileName"),
                exportJob.get("figure", "1"),
                float(exportJob.get("width", 469.47)),
                float(exportJob.get("scale", 1.0)),
                exportJob.get("heightRatio", False),
                exportJob.get("png", False),
            )

    def _fancyFigSize(self, scale: float, width: float, heightRatio: float = False) -> tuple[float, float]:
        """Create a fancy figure size compatible with matplotlib specs."""

        fig_width_pt = width
        inches_per_pt = 1.0 / 72.27  # Convert pt to inch
        golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio (you could change this)
        fig_width = fig_width_pt * inches_per_pt * scale  # width in inches
        fig_height = fig_width * (float(heightRatio) or golden_mean)  # height in inches
        fig_size = [fig_width, fig_height]
        return fig_size

    def exportFigure(
        self,
        fileName: str,
        figureID: int,
        width: float = 469.47,
        scale: float = 1.0,
        heightRatio: float = False,
        png: bool = False,
    ):
        """Export a figure.

        Parameters
        ----------
        fileName
            The filename.
        figureID
            The matplotlib figure ID.
        width
            The width in pt.
        scale
            Scale the figure.
        heightRatio.
            Set the height ratio or take golden mean if False.
        png
            Also export a .png figure.
        """

        fig, ax = self.figsWithAxes[figureID]
        fig.set_size_inches(self._fancyFigSize(scale, width, heightRatio))
        fig.tight_layout(pad=0.15)
        fig.savefig("{}.pgf".format(fileName))
        fig.savefig("{}.pdf".format(fileName))
        if png:
            fig.savefig("{}.png".format(fileName), dpi=400)

    def finalize(
        self,
    ):
        """Finalize and export plots."""

        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.configurePlotter()

            for fig, axes in self.figsWithAxes.values():
                fig.tight_layout()
                for axSpec, ax in axes.items():
                    leglabels = (
                        ax.get_legend_handles_labels()
                    )  # print legend only if labels were defined (ommit for contourplots)
                    if leglabels[0]:
                        ax.legend()
                    ax.relim()

        self.exportPlots()

    def show(
        self,
    ):
        """Show the plots!"""

        plt.show()
