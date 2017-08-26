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
from fe.utils.misc import stringDict
from fe.utils.meshtools import transferElsetResultsToElset, extractNodesFromElementSet
import fe.config.phenomena
import matplotlib.tri as mtri
from matplotlib import colors
import sympy as sp

documentation = {
        'figure': 'figure number, (default=1)',
        'axSpec': 'axis specification according to matplotlib syntax, (default=111)',
        'create=perNode': 'result per node is plotted in a meshplot',
        'create=perElement': 'result per element is plotted in a meshplot',
        'create=xyData': '2D Plot of results',
        }

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

class MeshPlot:

    def __init__(self,coordinates, elNodesIdxList, elCoordinatesList):
        self.coordinates = coordinates
        self.elCoordinatesList = elCoordinatesList
        self.xCoord = coordinates[:,0]
        self.yCoord = coordinates[:,1]
        self.xLimits = [-self.xCoord.max()*0.1, self.xCoord.max()*1.1]
        self.yLimits = [-self.yCoord.max()*0.1, self.yCoord.max()*1.1]
        self.TriangObj = Triangulation(self.xCoord, self.yCoord, elNodesIdxList)
        self.contourPlotScaling = 50
        self.userColorMap = 'coolwarm'
        
    def contourPlotFieldVariable(self, fieldValues, fig, ax, label):
        """ divide quad elements into two triangles and apply a constant field
        value for both triangles """
        resultPerTriElement = self.TriangObj.quadFieldToTriField(fieldValues)
        mapping = ax.tripcolor(self.xCoord, self.yCoord, self.TriangObj.triangleIdx, facecolors=resultPerTriElement, 
                                    cmap= self.userColorMap, norm=colors.Normalize(vmax=fieldValues.max(), vmin=fieldValues.min()) )
        cbar = fig.colorbar(mapping,fraction=0.046, pad=0.04)        
        cbar.set_label(label)
#        ax.set_xlim(self.xLimits)
#        ax.set_ylim(self.yLimits)

    def contourPlotNodalValues(self, z, fig, ax, label):
        """ divide quads into two triangles and apply a nodal value to the corner nodes """
        mapping = ax.tricontourf(self.TriangObj.triang, z, self.contourPlotScaling, cmap= self.userColorMap, norm=colors.Normalize(vmax=z.max(), vmin=z.min()) ) 
        cbar = fig.colorbar(mapping,fraction=0.046, pad=0.04)
        cbar.set_label(label)
        ax.set_xlim(self.xLimits)
        ax.set_ylim(self.yLimits)

    def plotNodeLabels(self, labels, ax):
        """ label nodes of elements """
        for label in labels:              
            ax.annotate('%i' % label, xy=self.coordinates[label-1,:], fontsize=6, textcoords='data')

    def plotMeshGrid(self, ax):
        """ plot grid of elements; so far only implemented for quads """
        for element in self.elCoordinatesList:
            ax.plot(np.append(element[:,0],element[0,0]),
                         np.append(element[:,1],element[0,1]), 'k', linewidth=0.3)
            
class OutputManager(OutputManagerBase):
    identification = "meshPlot"

    def __init__(self, name, definitionLines, jobInfo, modelInfo, fieldOutputController, journal, plotter):
        self.domainSize = jobInfo['domainSize']
        self.plotter = plotter
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
                    perNodeJob['fieldOutput'] = fieldOutputController.fieldOutputs[ definition['fieldOutput'] ]
                    if perNodeJob['fieldOutput'].type != 'perNode':
                        raise Exception('Meshplot: Please define perNode output on an nSet, not on a elSet!')
                    perNodeJob['label']  =          definition.get('label', '')
                    perNodeJob['axSpec'] =          int(definition.get('axSpec','111'))       
                    perNodeJob['figure'] =          int(definition.get('figure','1'))
                    perNodeJob['plotMeshGrid'] =    definition.get('plotMeshGrid', 'unDeformed')
                    if 'f(x)' in definition:
                        perNodeJob['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), definition['f(x)'] , 'numpy')
            
                    perNodeJob['plotNodeLabels'] =  definition.get('plotNodeLabels', False)
                    perNodeJob['dimensions'] =      fe.config.phenomena.getFieldSize(perNodeJob['fieldOutput'].field, self.domainSize)
                    self.perNodeJobs.append(perNodeJob)
                        
                elif varType == 'perElement':
                    perElementJob = {}
                    perElementJob['label']  =          definition.get('label', '')
                    perElementJob['axSpec'] =       int(definition.get('axSpec','111'))
                    perElementJob['figure'] =       int(definition.get('figure','1'))
                    perElementJob['fieldOutput'] =     fieldOutputController.fieldOutputs[ definition['fieldOutput'] ]
                    if 'f(x)' in definition:
                        perElementJob['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), definition['f(x)'] , 'numpy')
                            
                    perElementJob['plotMeshGrid'] = definition.get('plotMeshGrid', 'unDeformed')
                    perElementJob['plotNodeLabels'] =  definition.get('plotNodeLabels', False)
                    self.perElementJobs.append(perElementJob)
                    
                elif varType == 'xyData':
                    xyJob = {}
                    
                    if definition['x'] !=  'time':
                        xyJob['x'] = fieldOutputController.fieldOutputs[ definition['x'] ]
                    else:
                        xyJob['x'] = 'time'
                        
                    xyJob['y'] = fieldOutputController.fieldOutputs[ definition['y'] ]
                        
                    if 'f(x)' in definition:
                        xyJob['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), definition['f(x)'] , 'numpy')
                    if 'f(y)' in definition:
                        xyJob['f(y)'] = sp.lambdify ( sp.DeferredVector('y'), definition['f(y)'] , 'numpy')
                        
                    xyJob['figure'] =   int(definition.get('figure','1'))
                    xyJob['label'] =     definition.get('label', xyJob['y'].name)
                    xyJob['axSpec'] =   int(definition.get('axSpec','111'))
                    self.xyJobs.append(xyJob)
            
            if 'configFigure' in definition:
                self.configJobs.append(dict(definition))
            
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
        if self.perElementJobs or self.perNodeJobs:
            for element in self.elements.values():
                nodeArray = [node.coordinates for node in element.nodes][:4]
                nodeIdxArray = [nodeNumber.label-1 for nodeNumber in element.nodes[:]][:4]
                self.elCoordinatesList.append(np.asarray(nodeArray))
                self.elNodesIdxList.append(nodeIdxArray)
                
            self.meshPlot = MeshPlot(self.coordinateList, self.elNodesIdxList, self.elCoordinatesList)
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        pass
    
    def finalizeStep(self, U, P,):
        pass
        
    def finalizeJob(self, U, P,):
        for xyJob in self.xyJobs:
            
            y = xyJob['y'].getResultHistory()
            
            if xyJob['x'] == 'time':
                x = xyJob['y'].getTimeHistory()
            else:
                x = xyJob['x'].getResultHistory()
            
            if 'f(x)' in xyJob:
                x = xyJob['f(x)'] ( x )
            if 'f(y)' in xyJob:
                y = xyJob['f(y)'] ( y )
                
            self.plotter.plotXYData(x, y, xyJob['figure'], xyJob['axSpec'], xyJob )
            
        for perNodeJob in self.perNodeJobs:
            result = perNodeJob['fieldOutput'].getLastResult()
            
            fig = self.plotter.getFig(perNodeJob['figure'])
            ax =  self.plotter.getAx(perNodeJob['figure'] , perNodeJob['axSpec'])
            
            if 'f(x)' in perNodeJob:
                result = perNodeJob['f(x)'] (result)
    
            if perNodeJob['plotMeshGrid']=='unDeformed':
                self.meshPlot.plotMeshGrid( ax)
                
            if perNodeJob['plotNodeLabels']:
                self.meshPlot.plotNodeLabels(self.nodes.keys(),ax)
                
            result = np.squeeze(result)
            
            self.meshPlot.contourPlotNodalValues(result, fig, ax, perNodeJob['label'])

        for perElementJob in self.perElementJobs:

            fig = self.plotter.getFig(perElementJob['figure'])
            ax =  self.plotter.getAx(perElementJob['figure'] , perElementJob['axSpec'])
            self.meshPlot.plotMeshGrid(ax)
            
            if perElementJob['plotNodeLabels']:
                self.meshPlot.plotNodeLabels(self.nodes.keys(), ax)
                
            resultArray = perElementJob['fieldOutput'].getLastResult()
            
            if 'f(x)' in perElementJob:
                resultArray = perElementJob['f(x)'] (resultArray)
                
            if  perElementJob['fieldOutput'].elSetName != 'all':
                shape = ( len(self.elSets['all']), resultArray.shape[-1] ) if resultArray.ndim >= 2 else len(self.elSets['all'])
                resultsTarget = np.empty( shape  )
                resultsTarget[:] = np.nan
                transferElsetResultsToElset( self.elSets['all'], perElementJob['fieldOutput'].elSet, resultsTarget, resultArray )
                resultArray = resultsTarget
                
            self.meshPlot.contourPlotFieldVariable(resultArray, fig, ax, perElementJob['label'] )
            
        for configJob in self.configJobs:
            self.plotter.configAxes(**configJob)
        
        for saveJob in self.saveJobs:
            self.plotter.exportFigure(saveJob['fileName'], saveJob['figure'], saveJob['width'], saveJob['scale'], saveJob['heightRatio'], saveJob['png'])
       
