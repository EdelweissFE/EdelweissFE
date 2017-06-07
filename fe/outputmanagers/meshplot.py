#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:56:11 2017

@author: c8441146
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
import sympy as sp

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
    def __init__(self, coordinates, elNodesIdxList, elCoordinatesList):
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
        
        self.coordinates = coordinates
        self.elCoordinatesList = elCoordinatesList
        self.xCoord = coordinates[:,0]
        self.yCoord = coordinates[:,1]
        self.xLimits = [-self.xCoord.max()*0.1, self.xCoord.max()*1.1]
        self.yLimits = [-self.yCoord.max()*0.1, self.yCoord.max()*1.1]
        self.TriangObj = Triangulation(self.xCoord, self.yCoord, elNodesIdxList)
    
    def setFigAxesLabel(self, axSpec=111, figureID=0, userLabel='', cmap='coolwarm'):
        """ create a figure with axes if it doesn't exist so far """
        if not figureID in self.figsWithAxes:
            self.figsWithAxes[figureID] = (plt.figure(figureID), {})
        self.fig, axes = self.figsWithAxes[figureID] 
        self.ax = self.fig.add_subplot(axSpec)
        self.userLabel = userLabel
        self.userColorMap = cmap

    def contourPlotFieldVariable(self, fieldValues):
        """ divide quad elements into two triangles and apply a constant field
        value for both triangles """
        resultPerTriElement = self.TriangObj.quadFieldToTriField(fieldValues)
        mapping = self.ax.tripcolor(self.xCoord, self.yCoord, self.TriangObj.triangleIdx, facecolors=resultPerTriElement, 
                                    cmap= self.userColorMap, norm=colors.Normalize(vmax=fieldValues.max(), vmin=fieldValues.min()) )
        cbar = self.fig.colorbar(mapping)
        cbar.set_label(self.userLabel)
        self.ax.set_xlim(self.xLimits)
        self.ax.set_ylim(self.yLimits)

        
    def contourPlotNodalValues(self, z):
        """ divide quads into two triangles and apply a nodal value to the corner nodes """
        mapping = self.ax.tricontourf(self.TriangObj.triang, z, self.contourPlotScaling, cmap= self.userColorMap, norm=colors.Normalize(vmax=z.max(), vmin=z.min()) ) 
        cbar = self.fig.colorbar(mapping)
        cbar.set_label(self.userLabel)
        self.ax.set_xlim(self.xLimits)
        self.ax.set_ylim(self.yLimits)
        
    def plotNodeLabels(self, labels):
        """ label nodes of elements """
        for label in labels:              
            self.ax.annotate('%i' % label, xy=self.coordinates[label-1,:], fontsize=6, textcoords='data')
            
    def plotMeshGrid(self):
        """ plot grid of elements; so far only implemented for quads """
        for element in self.elCoordinatesList:
            self.ax.plot(np.append(element[:,0],element[0,0]),
                         np.append(element[:,1],element[0,1]), 'k', linewidth=0.3)

    def fancyFigSize(self, scale, width, heightRatio=False):
        fig_width_pt = width                       
        inches_per_pt = 1.0/72.27                       # Convert pt to inch
        golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
        fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
        fig_height = fig_width* (float(heightRatio) or golden_mean)              # height in inches
        fig_size = [fig_width,fig_height]
        return fig_size
    
    def exportFigure(self, ):
        pass
    
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
        
        # write List of nodeLabels
        self.labelList = np.asarray([nodeNumber for nodeNumber in self.nodes.keys()])
        # write List of node coordiantes
        self.coordinateList = np.asarray([node.coordinates for node in self.nodes.values()])    
        # write list of element coordinates with 4x2 arrays (xCol, yCol)
        self.elCoordinatesList = []
        # write list of node indices for each element relevant for the meshplot output
        # in case of an 8-node element only the 4 first nodes are relevant
        self.elNodesIdxList = []
        for element in self.elements.values():
            nodeArray = [node.coordinates for node in element.nodes][:4]
            nodeIdxArray = [nodeNumber.label-1 for nodeNumber in element.nodes[:]][:4]
            self.elCoordinatesList.append(np.asarray(nodeArray))
            self.elNodesIdxList.append(nodeIdxArray)
        
        self.perNodeJobs = []
        self.perElementJobs = []
        
        for defLine in definitionLines:
            definition = stringDict(defLine)
              
            if 'create' in definition:
                varType = definition['create']
                
                if varType == 'perNode':
                    perNodeJob = {}
                    perNodeJob['name'] =            definition.get('name', 'defaultJob')
                    field =                         definition['field']
                    perNodeJob['result'] =          definition.get('result')
                    perNodeJob['axSpec'] =          int(definition.get('axSpec','111'))       
                    perNodeJob['figure'] =          int(definition.get('figure','1'))
                    perNodeJob['plotMeshGrid'] =    definition.get('plotMeshGrid', 'unDeformed')
                    perNodeJob['indices'] =         np.asarray([node.fields[field] for node in self.nodes.values()]).ravel()
                    
                    f = definition.get('f(x)', 'x')
                    perNodeJob['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), f , 'numpy')
            
                    perNodeJob['resultIndices'] =   [node.fields[definition['field']] for node in self.nodes.values()]
                    perNodeJob['plotNodeLabels'] =  definition.get('plotNodeLabels', False)
                    perNodeJob['dimensions'] =      fe.config.phenomena.getFieldSize(field, self.domainSize)
                    self.perNodeJobs.append(perNodeJob)
                        
                if varType == 'perElement':
                    perElementJob = {}
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
        
        # Create an Instance of Plotter class
        self.plotObj = Plotter(self.coordinateList, self.elNodesIdxList, self.elCoordinatesList)
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass

    def finalizeStep(self, U, P,):
        pass
        
    def finalizeJob(self, U, P,):
        for perNodeJob in self.perNodeJobs:
            self.plotObj.setFigAxesLabel(perNodeJob['axSpec'], perNodeJob['figure'], perNodeJob['name'])
            location = U if perNodeJob['result'] == 'U' else P
                
            indices = perNodeJob['resultIndices']
            result = np.asarray([perNodeJob['f(x)'](row) for row in location[indices]])
            if perNodeJob['plotMeshGrid']=='unDeformed':
                self.plotObj.plotMeshGrid()
            if perNodeJob['plotNodeLabels']:
                self.plotObj.plotNodeLabels(self.nodes.keys())
            self.plotObj.contourPlotNodalValues(result)

        for perElementJob in self.perElementJobs:
            self.plotObj.setFigAxesLabel(perElementJob['axSpec'], perElementJob['figure'], perElementJob['name'])
            self.plotObj.plotMeshGrid()
            if perElementJob['plotNodeLabels']:
                self.plotObj.plotNodeLabels(self.nodes.keys())
            
            resultArray = np.asarray([el.getPermanentResultPtr(**perElementJob) for el in  self.elements.values()])
            resultIndex = perElementJob['idxStart']
            if perElementJob['result'] == 'sdv':
                self.plotObj.contourPlotFieldVariable(resultArray)
            else:
                self.plotObj.contourPlotFieldVariable(np.asarray([ result[resultIndex] for result in resultArray  ]))
              
        self.plotObj.show()
        
    def finalizeIncrement(self, U, P, increment):
        pass
                
