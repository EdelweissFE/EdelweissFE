#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Magdalena

Module meshplot divided into classes:
    * Triangulation:
        creates triangles out of rectangles in order to use matplotlib plotting
        options
    * Plotter:
        creates figures, axes, grid, labels
    * Outputmanager:
        creats the plotting specific for the defined keyword lines
"""
from fe.outputmanagers.outputmanagerbase import OutputManagerBase
import numpy as np
import os 
import sys
from fe.utils.misc import stringDict
import fe.config.phenomena
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import colors
import matplotlib
import itertools
import sympy as sp

documentation = {
        'figure': 'figure number, (default=1)',
        'axSpec': 'axis specification according to matplotlib syntax, (default=111)',
        'create=perNode': 'result per node is plotted in a meshplot',
        'create=perElement': 'result per element is plotted in a meshplot',
        'create=xyData': '2D Plot of results',
        }

defaultMarkerCycle = itertools.cycle(('o', 'v', 'D', 's', '^'))
defaultLinePlotColorCycle = itertools.cycle(('k'))
defaultScatterPlotColorCycle = itertools.cycle(('k', 'r', 'b', 'g'))
defaultLineStyleCycle = itertools.cycle(('-', '0-3-2', '0-3-1-1-1', '0-1-1'))

class Triangulation:
    """ class that provides the division of quadrilateral elements into triangles """
    def __init__(self, xCoord, yCoord, elNodesIdxList):
        self.triangleIdx, self.triang = self.quadIdxToTriIdx(xCoord, yCoord, elNodesIdxList)
        
    def quadIdxToTriIdx(self, xCoord, yCoord, elementIdxMatrix):
        triangleIdx = np.asarray([ [list(elIdx[:3]) , [elIdx[2],elIdx[3],elIdx[0] ]  ] for elIdx in elementIdxMatrix  ]).reshape(-1,3)      
        triang = mtri.Triangulation(xCoord,yCoord,triangleIdx)
        return triangleIdx, triang 
        
    def quadFieldToTriField(self, fieldValues):
        return np.asarray([2*[fieldValues]]).reshape(-1,1, order='f').flatten()

class Plotter:
    def __init__(self):
        self.rcParams = {
                "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
                "text.usetex": True,                # use LaTeX to write all text
                'text.latex.preamble':[r"\usepackage{mathpazo}",
                                       r"\usepackage{siunitx}",
                                        ],
                "font.family": "serif",
                "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
                "font.sans-serif": [],
                "font.monospace": [],
                "axes.labelsize": 10,               # LaTeX default is 10pt font.
                "font.size": 10,
                "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
                "legend.numpoints":1,
                'legend.labelspacing':.2,
                "xtick.labelsize": 8,
                "ytick.labelsize": 8,
                'lines.linewidth': 1,
                'figure.dpi':200,
                'lines.markeredgewidth': 0.4,
                'lines.markersize':4
                }
        plt.close('all')
        
        self.contourPlotScaling = 50
        
        if(os.path.isfile('plotterRcParams.py')):  
            print("found plotterRcParams.py")
            sys.path.append(os.getcwd())
            import plotterRcParams
            
            self.rcParams.update(plotterRcParams.rcParams)
#        
        self.figsWithAxes = {} # {figID : (figure, {axesDict})}
        matplotlib.rcParams.update(self.rcParams)
        
    def initMeshPlot(self,coordinates, elNodesIdxList, elCoordinatesList):
        self.coordinates = coordinates
        self.elCoordinatesList = elCoordinatesList
        self.xCoord = coordinates[:,0]
        self.yCoord = coordinates[:,1]
        self.xLimits = [-self.xCoord.max()*0.1, self.xCoord.max()*1.1]
        self.yLimits = [-self.yCoord.max()*0.1, self.yCoord.max()*1.1]
        self.TriangObj = Triangulation(self.xCoord, self.yCoord, elNodesIdxList)
    
    def setFigLabelAxes(self, figureID=0, axSpec=111, userLabel='', cmap='coolwarm'):
        """ create a figure with label if it doesn't exist so far """
        if not figureID in self.figsWithAxes:
            self.figsWithAxes[figureID] = (plt.figure(figureID),{})
        
        self.fig, self.axes = self.figsWithAxes[figureID]
        if not axSpec in self.axes:
            self.axes[axSpec] = axSpec
            self.figsWithAxes[figureID][1][axSpec] =  self.fig.add_subplot(axSpec)
        self.userLabel = userLabel
        self.userColorMap = cmap

    def contourPlotFieldVariable(self, fieldValues, axSpec=111):
        """ divide quad elements into two triangles and apply a constant field
        value for both triangles """
        resultPerTriElement = self.TriangObj.quadFieldToTriField(fieldValues)
        mapping = self.axes[axSpec].tripcolor(self.xCoord, self.yCoord, self.TriangObj.triangleIdx, facecolors=resultPerTriElement, 
                                    cmap= self.userColorMap, norm=colors.Normalize(vmax=fieldValues.max(), vmin=fieldValues.min()) )
        cbar = self.fig.colorbar(mapping)        
        cbar.set_label(self.userLabel)
        self.axes[axSpec].set_xlim(self.xLimits)
        self.axes[axSpec].set_ylim(self.yLimits)

    def contourPlotNodalValues(self, z, axSpec=111):
        """ divide quads into two triangles and apply a nodal value to the corner nodes """
        mapping = self.axes[axSpec].tricontourf(self.TriangObj.triang, z, self.contourPlotScaling, cmap= self.userColorMap, norm=colors.Normalize(vmax=z.max(), vmin=z.min()) ) 
        cbar = self.fig.colorbar(mapping)
        cbar.set_label(self.userLabel)
        self.axes[axSpec].set_xlim(self.xLimits)
        self.axes[axSpec].set_ylim(self.yLimits)
        
    def plotXYData(self, x, y, figureID=1, axSpec=111, plotOptions=None):
        """ plots a single curve """        
        plotDefinition = {}
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
        self.figsWithAxes[figureID][1][axSpec].plot(x,y, **plotDefinition)   
        
    def plotNodeLabels(self, labels, axSpec=111):
        """ label nodes of elements """
        for label in labels:              
            self.axes[axSpec].annotate('%i' % label, xy=self.coordinates[label-1,:], fontsize=6, textcoords='data')
            
    def plotMeshGrid(self, axSpec=111):
        """ plot grid of elements; so far only implemented for quads """
        for element in self.elCoordinatesList:
            self.axes[axSpec].plot(np.append(element[:,0],element[0,0]),
                         np.append(element[:,1],element[0,1]), 'k', linewidth=0.3)

        
    def configAxes(self, configEntry):
        ax = self.figsWithAxes[configEntry['figure']][1][configEntry['axSpec']]
        if configEntry["xLimits"] is not None: 
            limits = [float(x) for x in configEntry["xLimits"].split('_')] 
            ax.set_xlim(limits)
        if configEntry["yLimits"] is not None: 
            limits = [float(x) for x in configEntry["yLimits"].split('_')] 
            ax.set_ylim(limits)
        ax.set_xlabel(configEntry['xLabel'])
        ax.set_ylabel(configEntry['yLabel'])
        if configEntry.get('flipX', 'False') == 'True':   ax.invert_xaxis()
        if configEntry.get('flipY', 'False') == 'True':   ax.invert_yaxis()

    def fancyFigSize(self, scale, width, heightRatio=False):
        fig_width_pt = width                       
        inches_per_pt = 1.0/72.27                       # Convert pt to inch
        golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
        fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
        fig_height = fig_width* (float(heightRatio) or golden_mean)              # height in inches
        fig_size = [fig_width,fig_height]
        return fig_size
    
    def exportFigure(self, fileName, figure, width=469.47, scale=1.0, heightRatio=False, png=False):
        fig, ax = self.figsWithAxes[figure] 
        fig.set_size_inches(self.fancyFigSize(scale, width, heightRatio))
        fig.tight_layout(pad=0.15)
        fig.savefig('{}.pgf'.format(fileName))
        fig.savefig('{}.pdf'.format(fileName))
        if png:
            fig.savefig('{}.png'.format(fileName), dpi=400)
    
    def show(self,):
        plt.show()
    
            
class OutputManager(OutputManagerBase):
    identification = "meshPlot"
    #printTemplate = "node {:}, {:} {:} {:}: {:}"
    resultFunctions = {'1' : lambda x: x[0],
                       '2' : lambda x: x[1],
                       '3' : lambda x: x[2],
                       'all' : lambda x: x,
                       'sum' : lambda x: np.sum(x),
                       'mean' : lambda x: np.mean(x),
                       'magnitude' : lambda x: np.linalg.norm(x),
                       }

    def __init__(self, name, definitionLines, jobInfo, modelInfo, journal):
        self.finishedSteps =0
        self.domainSize = jobInfo['domainSize']
        self.journal = journal

        self.nodes = modelInfo['nodes']
        self.elements = modelInfo['elements']
        self.elSets = modelInfo['elementSets']
        self.nSets = modelInfo['nodeSets']
        
        # write List of nodeLabels
        self.labelList = np.asarray([nodeNumber for nodeNumber in self.nodes.keys()])
        # write List of node coordiantes
        self.coordinateList = np.asarray([node.coordinates for node in self.nodes.values()])    
        # write list of element coordinates with 4x2 arrays (xCol, yCol)
        self.elCoordinatesList = []
        # write list of node indices for each element relevant for the meshplot output
        # in case of an 8-node element only the 4 first nodes are relevant
        self.elNodesIdxList = []
        
        self.perNodeJobs = []
        self.perElementJobs = []
        self.configJobs = []
        self.xyJobs = []
        self.saveJobs = []
        
        for defLine in definitionLines:
            definition = stringDict(defLine)
            
            if 'create' in definition:
                varType = definition['create']
                
                if varType == 'perNode':
                    perNodeJob = {}
                    perNodeJob['name'] =            definition.get('name', 'defaultJob')
                    perNodeJob['field'] =           definition['field']
                    perNodeJob['result'] =          definition.get('result')
                    perNodeJob['axSpec'] =          int(definition.get('axSpec','111'))       
                    perNodeJob['figure'] =          int(definition.get('figure','1'))
                    perNodeJob['plotMeshGrid'] =    definition.get('plotMeshGrid', 'unDeformed')
                    
                    f = definition.get('f(x)', 'x')
                    perNodeJob['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), f , 'numpy')
            
                    perNodeJob['resultIndices'] =   [node.fields[definition['field']] for node in self.nodes.values()]
                    perNodeJob['plotNodeLabels'] =  definition.get('plotNodeLabels', False)
                    perNodeJob['dimensions'] =      fe.config.phenomena.getFieldSize(perNodeJob['field'], self.domainSize)
                    self.perNodeJobs.append(perNodeJob)
                        
                elif varType == 'perElement':
                    perElementJob = {}
                    perElementJob['elSet']  =       definition.get('elSet','all')
                    perElementJob['axSpec'] =       int(definition.get('axSpec','111'))
                    perElementJob['figure'] =       int(definition.get('figure','1'))
                    perElementJob['result'] =       definition['result']
                    if 'index' in definition:
                        idcs = definition['index']
                        if ':' in definition['index']:
                            idcs=[int (i) for i in idcs.split(':')]
                            perElementJob['idxStart'], perElementJob['idxStop'] = idcs
                        else:
                            idx = int (idcs )
                            perElementJob['idxStart'], perElementJob['idxStop'] = idx, idx+1
                    perElementJob['plotMeshGrid'] = definition.get('plotMeshGrid', 'unDeformed')
                    perElementJob['plotNodeLabels'] =  definition.get('plotNodeLabels', False)
                    perElementJob['gaussPt'] =      int(definition['gaussPt'])
                    perElementJob['name'] =         definition.get('name', perElementJob['result'])
                    self.perElementJobs.append(perElementJob)
                    
                elif varType == 'xyData':
                    xyJob = {}
                    for key, value in definition.items():
                        xyJob[key] = value
                    
                    xyJob['figure'] = int(definition.get('figure','1'))
                    xyJob['name'] =         definition.get('name', '')
                    xyJob['axSpec'] = int(definition.get('axSpec','111'))
                    xyJob['xResult'] = definition.get('xResult')
                    xyJob['area'] = definition.get('area', False)
                    if xyJob['xResult'][0] == 'U' or 'P': 
                        xyJob['xField'] = definition.get('xField')
                        xyJob['xDoF'] = int(xyJob['xResult'][1])
                        xyJob['xnSet'] = definition.get('xnSet') or definition.get('nSet') 
                        f = definition.get('f(x)', 'x')
                        xyJob['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), f , 'numpy')
                                    #Initialize result array
                        xyJob['xResultIndices'] = np.asarray([node.fields[xyJob['xField']] for node in self.nSets[xyJob['xnSet']]])
                        
                    xyJob['yResult'] = definition.get('yResult')
                    if xyJob['yResult'][0] == 'U' or 'P': 
                        xyJob['yField'] = definition.get('yField')
                        xyJob['yDoF'] = int(xyJob['yResult'][1])
                        xyJob['ynSet'] = definition.get('ynSet') or definition.get('nSet') 
                        f = definition.get('f(y)', 'y')
                        xyJob['f(y)'] = sp.lambdify ( sp.DeferredVector('y'), f , 'numpy')
                        xyJob['yResultIndices'] = np.asarray([node.fields[xyJob['yField']] for node in self.nSets[xyJob['ynSet']]])   
                    xyJob['xList'] = []
                    xyJob['yList'] = []
                    self.xyJobs.append(xyJob)
            
            if 'configFigure' in definition:
                configJob = {}
                configJob['figure'] = int(definition.get('figure','1'))
                configJob['axSpec'] = int(definition.get('axSpec','111'))
                configJob['xLimits'] = definition.get('xLimits',None)
                configJob['yLimits'] = definition.get('yLimits',None)
                configJob['xLabel'] = definition.get('xLabel',None)
                configJob['yLabel'] = definition.get('yLabel',None)
                configJob['flipX'] = definition.get('flipX',False)
                configJob['flipY'] = definition.get('flipY',False)
                configJob['marker'] = definition.get('marker',None)
                self.configJobs.append(configJob)
            
            if 'saveFigure' in definition:
                saveJob = {}
                saveJob['figure'] = int(definition.get('figure','1'))
                saveJob['fileName'] = definition.get('fileName')
                saveJob['width'] = definition.get('width',469.47)
                saveJob['scale'] = definition.get('scale', 1.0)
                saveJob['heightRatio'] = definition.get('heightRatio', False)
                saveJob['png'] = definition.get('png', False)
                self.saveJobs.append(saveJob)        

        # Initialize instance of plotterclass                 
        self.plotObj = Plotter()
        if self.perElementJobs or self.perNodeJobs:
            for element in self.elements.values():
                nodeArray = [node.coordinates for node in element.nodes][:4]
                nodeIdxArray = [nodeNumber.label-1 for nodeNumber in element.nodes[:]][:4]
                self.elCoordinatesList.append(np.asarray(nodeArray))
                self.elNodesIdxList.append(nodeIdxArray)
            self.plotObj.initMeshPlot(self.coordinateList, self.elNodesIdxList, self.elCoordinatesList)
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        for xyJob in self.xyJobs:
            xLocation = U if xyJob['xResult'][0] == 'U' else P
            yLocation = U if xyJob['yResult'][0] == 'U' else P
            
            # indices dependent on field variable
            xIndices = xyJob['xResultIndices']
            yIndices = xyJob['yResultIndices']
            # indices dependent on function f(x)
            
            xResult =xyJob['f(x)'](xLocation[xIndices][:,xyJob['xDoF']-1])
            yResult = xyJob['f(y)'](yLocation[yIndices][:,xyJob['yDoF']-1])
            xyJob['xList'].append(xResult)
            xyJob['yList'].append(yResult)

    def finalizeStep(self, U, P,):
        for xyJob in self.xyJobs:
            xLocation = U if xyJob['xResult'][0] == 'U' else P
            yLocation = U if xyJob['yResult'][0] == 'U' else P
            # indices dependent on field variable
            xIndices = xyJob['xResultIndices']
            yIndices = xyJob['yResultIndices']
            # indices dependent on function f(x)
            xResult =xyJob['f(x)'](xLocation[xIndices][:,xyJob['xDoF']-1])
            yResult = xyJob['f(y)'](yLocation[yIndices][:,xyJob['yDoF']-1])
            
            xyJob['xList'].append(xResult)
            xyJob['yList'].append(yResult)
        
    def finalizeJob(self, U, P,):
        for xyJob in self.xyJobs:
            self.plotObj.setFigLabelAxes(xyJob['figure'], xyJob['axSpec'], xyJob['name'])
            self.plotObj.plotXYData(xyJob['xList'], xyJob['yList'], xyJob['figure'], xyJob['axSpec'], xyJob )
            if xyJob['area']:
                xyJobArea = np.trapz( xyJob['yList'], x=xyJob['xList'])
                self.plotObj.axes[xyJob['axSpec']].annotate('A='+str(xyJobArea),xy=(np.mean(xyJob['xList']),np.mean(xyJob['yList'])))
                self.plotObj.axes[xyJob['axSpec']].fill(xyJob['xList'], xyJob['yList'], color='gray', alpha=0.5 )
            
        for perNodeJob in self.perNodeJobs:
            self.plotObj.setFigLabelAxes(perNodeJob['figure'],perNodeJob['axSpec'], perNodeJob['name'])
            location = U if perNodeJob['result'] == 'U' else P
            # indices dependent on field variable
            indices = perNodeJob['resultIndices']
            result = []
            for row in indices:
                res = []
                for idx in row:
                    res.append(location[idx])
                result.append(perNodeJob['f(x)'](np.asarray(res)))
            
            result = np.asarray(result)
#            result = np.asarray([perNodeJob['f(x)'](row) for row in location[indices]])  # DID NOT WORK OUT
            if perNodeJob['plotMeshGrid']=='unDeformed':
                self.plotObj.plotMeshGrid(perNodeJob['axSpec'],)
            if perNodeJob['plotNodeLabels']:
                self.plotObj.plotNodeLabels(self.nodes.keys(),perNodeJob['axSpec'])
            self.plotObj.contourPlotNodalValues(result,perNodeJob['axSpec'])

        for perElementJob in self.perElementJobs:
            self.plotObj.setFigLabelAxes(perElementJob['figure'], perElementJob['axSpec'], perElementJob['name'])
            self.plotObj.plotMeshGrid(perElementJob['axSpec'])
            if perElementJob['plotNodeLabels']:
                self.plotObj.plotNodeLabels(self.nodes.keys(), perElementJob['axSpec'])
            #Initialize result array
            resultArray = np.asarray([np.nan for el in  self.elements.values()])
            #Fill only entries where element results are present
            for el in self.elSets[perElementJob['elSet']]:
                resultArray[el.elNumber-1] = el.getPermanentResultPtr(**perElementJob)

            resultIndex = perElementJob['idxStart']
            if perElementJob['result'] == 'sdv':
                self.plotObj.contourPlotFieldVariable(resultArray, perElementJob['axSpec'])
            else:
                self.plotObj.contourPlotFieldVariable(np.asarray([ result[resultIndex] for result in resultArray  ]), perElementJob['axSpec'])
        
        for configJob in self.configJobs:
            self.plotObj.configAxes(configJob)
        
        for saveJob in self.saveJobs:
            self.plotObj.exportFigure(saveJob['fileName'], saveJob['figure'], saveJob['width'], saveJob['scale'], saveJob['heightRatio'], saveJob['png'])
        self.plotObj.show()
       
