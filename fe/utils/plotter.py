#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:30:36 2017

@author: matthias
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import itertools
from fe.utils.misc import stringDict

defaultMarkerCycle = itertools.cycle(('o', 'v', 'D', 's', '^'))
defaultLinePlotColorCycle = itertools.cycle(('k'))
defaultScatterPlotColorCycle = itertools.cycle(('k', 'r', 'b', 'g'))
defaultLineStyleCycle = itertools.cycle(('-', '0-3-2', '0-3-1-1-1', '0-1-1'))

class Plotter:
    """
    A Unified Plotter, which can be accessed and used by all outputmanagers
    """
    def __init__(self, journal, inputfile):
        self.journal = journal
        self.rcParams = {
                "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
                "text.usetex": True,                # use LaTeX to write all text
                'text.latex.preamble':r"\usepackage{mathpazo} \usepackage{siunitx}",
                "font.family": "serif",
#                "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
#                "font.sans-serif": [],
                "font.monospace": [],
                "axes.labelsize": 10,               # LaTeX default is 10pt font.
                "font.size": 10,
                "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
                "legend.numpoints":1,
                'legend.labelspacing':.2,
                "xtick.labelsize": 8,
                "ytick.labelsize": 8,
                'lines.linewidth': 1,
                'lines.markeredgewidth': 0.4,
                'lines.markersize':4
                }
        plt.close('all')
        
        
        if(os.path.isfile('plotterRcParams.py')):  
            print("found plotterRcParams.py")
            sys.path.append(os.getcwd())
            import plotterRcParams
            self.rcParams.update(plotterRcParams.rcParams)
            
        self.figsWithAxes = {} # {figID : (figure, {axesDict})}
        matplotlib.rcParams.update(self.rcParams)
        
        self.configurationLines = [ stringDict(c) for configEntry in  inputfile['*configurePlots'] for c in configEntry['data']   ]
        self.exportJobs = [ stringDict(c) for configEntry in  inputfile['*exportPlots'] for c in configEntry['data']   ]
        
    def getAx(self, figureID=0, axSpec=111):
        """ create a figure with ax if it doesn't exist so far """
        if not figureID in self.figsWithAxes:
            self.figsWithAxes[figureID] = (plt.figure(figureID),{})
        
        fig, axes = self.figsWithAxes[figureID]
        axSpec = axSpec
        if not axSpec in axes:
            self.figsWithAxes[figureID][1][axSpec] =  fig.add_subplot(int(axSpec))
            self.figsWithAxes[figureID][1][axSpec].grid(True)
        
        return self.figsWithAxes[figureID][1][axSpec]
    
    def getFig(self, figureID=0):
        """ create a figure doesn't exist so far """
        if not figureID in self.figsWithAxes:
            self.figsWithAxes[figureID] = (plt.figure(figureID),{})
        
        return self.figsWithAxes[figureID][0]
    
        
    def plotXYData(self, x, y, figureID=1, axSpec=111, plotOptions=None):
        """ plots a single curve """        
        
        ax = self.getAx(figureID, axSpec)
        
        plotDefinition = {}
        if 'label' in plotOptions: plotDefinition['label'] = plotOptions['label']
        if plotOptions.get('plotType', 'linePlot') == "scatter":
            plotDefinition['marker'] =          plotOptions.get('marker', None) or next(defaultMarkerCycle)
            if 'ms' in plotOptions:
                plotDefinition['markersize'] =      plotOptions.get('ms', None)
            plotDefinition['markerfacecolor'] = plotOptions.get('markerfacecolor', 'None') 
            plotDefinition['markeredgecolor'] = plotOptions.get('c', False) or next(defaultScatterPlotColorCycle)
            plotDefinition['linestyle'] = ''
        else:
            plotDefinition['c'] = plotOptions.get('c', False) or next(defaultLinePlotColorCycle) 
            plotDefinition['linestyle'] = plotOptions.get('ls', False) or next(defaultLineStyleCycle)
            lsParts = plotDefinition['linestyle'].split('-')
            if len(lsParts) >= 3:
                plotDefinition['linestyle'] = (float(lsParts[0]), [ int(onOff) for onOff in lsParts[1:] ])
        ax.plot(x,y, **plotDefinition)   
        
    def configurePlotter(self):
        
        for configEntry in self.configurationLines:
        
            ax = self.figsWithAxes[ ( configEntry['figure'] ) ][1][ configEntry['axSpec'] ]
            
            if "xLimits" in configEntry: 
                limits = [float(x) for x in configEntry["xLimits"].split('_')] 
                ax.set_xlim(limits)
                
            if "yLimits" in configEntry: 
                limits = [float(x) for x in configEntry["yLimits"].split('_')] 
                ax.set_ylim(limits)
                
            if 'xLabel' in configEntry: ax.set_xlabel(configEntry['xLabel'])
            if 'yLabel' in configEntry: ax.set_ylabel(configEntry['yLabel'])
            if 'flipX' in configEntry: ax.invert_xaxis()
            if 'flipY' in configEntry: ax.invert_yaxis()
            if 'aspect' in configEntry: ax.set_aspect(configEntry['aspect'])
            if 'grid' in configEntry: ax.grid()

    def exportPlots(self):
        for exportJob in self.exportJobs:
            self.exportFigure(exportJob.get('fileName'),
                          exportJob.get('figure','1'), 
                          float( exportJob.get('width',469.47) ), 
                          float( exportJob.get('scale', 1.0) ), 
                          exportJob.get('heightRatio', False), 
                          exportJob.get('png', False))
            
    def fancyFigSize(self, scale, width, heightRatio=False):
        fig_width_pt = width                       
        inches_per_pt = 1.0/72.27                       # Convert pt to inch
        golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
        fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
        fig_height = fig_width* (float(heightRatio) or golden_mean)              # height in inches
        fig_size = [fig_width,fig_height]
        return fig_size
    
    def exportFigure(self, fileName, figureID, width=469.47, scale=1.0, heightRatio=False, png=False):
        fig, ax = self.figsWithAxes[figureID]
        fig.set_size_inches(self.fancyFigSize(scale, width, heightRatio))
        fig.tight_layout(pad=0.15)
        fig.savefig('{}.pgf'.format(fileName))
        fig.savefig('{}.pdf'.format(fileName))
        if png:
            fig.savefig('{}.png'.format(fileName), dpi=400)
    
    def finalize(self,):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.configurePlotter()
            
            for fig, axes in self.figsWithAxes.values():
                fig.tight_layout()
                for axSpec, ax in axes.items():
                    leglabels = ax.get_legend_handles_labels() # print legend only if labels were defined (ommit for contourplots)
                    if leglabels[0]:
                        ax.legend()
                    ax.relim()
                    
        self.exportPlots()

    def show(self,):
        plt.show()
    
