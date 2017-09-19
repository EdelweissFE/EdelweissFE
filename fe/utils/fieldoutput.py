#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 14:57:48 2017

@author: matthias

FieldOutputs store all kind of analysis results, 
and are defined via the keyword *fieldOutput.
All fieldOutputs are accessable to all outputmanagers at the end of each 
increment, step and job.
Furthermore, they can be exported to *.csv files at the end of the analysis job.

ATTENTION: 
    If the results are exported to a .csv file with enabled "saveHistory", 
    the time History is automatically appended to the .csv file"

Datalines:
"""

documentation={'name':'name of the fieldOutput',
               'nSet|elSet|node|element': 'entity, for which the fieldOutput is defined',
               'result' : 'e.g., U, P, stress, strain ...',
               'gaussPt': 'for element based fieldOutputs only, integers or slices',
               'f(x)': '(optional), apply math (in each increment)',
               'saveHistory': '(optional), save complete History or only last (increment) result. Default: True (node, element) and False (nSet, elSet)',
               'export':'(optional), export the fieldOutput to a file at the end of the job',
               'f_export(x)': '(optional), apply math on the final result (table)'}

import numpy as np
from fe.utils.misc import stringDict, strToRange, isInteger
from fe.utils.meshtools import  extractNodesFromElementSet
from fe.utils.math import createMathExpression
from fe.utils.elementresultcollector import ElementResultCollector
from time import time as getCurrentTime

class FieldOutput:
    """
    Entity of a fieldOutput request
    """
    def __init__(self, modelInfo, definition, journal):
        self.name = definition['name']
        self.journal = journal
        self.timeTotal = 0.0
        # determination of type:
            # perNode, perElement, perNodeSet, perElset...
         
        if 'nSet' in definition:
            self.type = 'perNodeSet'
            self.nSet = modelInfo['nodeSets'] [ definition['nSet'] ]
            self.nSetName = definition['nSet']
            
        elif 'elSet' in definition:
            if  definition['result'] == 'U' or  definition['result'] == 'P':
                self.journal.message('Converting elSet {:} to a nSet due to requested nodal results'.format(definition['elSet']), self.name )
                # an elset was given, but in fact a nodeset was 'meant': we extract the nodes of the elementset!
                self.type = 'perNodeSet'
                self.nSet = extractNodesFromElementSet( modelInfo['elementSets'] [ definition['elSet'] ] )
                self.nSetName = definition['elSet']
            else:
                # it's really an elSet job
                self.type= 'perElementSet'
                self.elSet = modelInfo['elementSets'] [ definition['elSet'] ]
                self.elSetName =  definition['elSet']
            
        elif 'node' in definition:
            self.type = 'perNode'
            self.node = modelInfo['nodes'] [ int(definition['node']) ]
            
        elif 'element' in definition:
            self.type  = 'perElement'
            self.elSet = [ modelInfo['elements'] [ int(definition['element']) ] ]
        else:
            raise Exception('invalid field output requested: ' + definition['name'] )
            
        # save history, or only current value(s) ?
        if self.type == 'perNode' or self.type == 'perElement':
            self.appendResults = definition.get('saveHistory', True)
        else:
            self.appendResults = definition.get('saveHistory', False)
            
        self.result = [] if self.appendResults else None
            
        if self.type == 'perNode' or self.type == 'perNodeSet':
            self.resultVector =     definition['result']
            self.field =            definition['field']
            
            if self.type == 'perNodeSet':
                self.resultIndices =  np.array([ [ n.fields[ self.field ] ] for n in self.nSet], dtype=np.int).ravel()
                
        elif self.type == 'perElement' or self.type == 'perElementSet':
            
#            requestDictForElement = {}
#            requestDictForElement['result'] = definition['result']
            
            gpt = definition['gaussPt']
            gaussPts = strToRange(gpt) if not isInteger(gpt) else [ int (gpt)  ]
            
            self.elementResultCollector = ElementResultCollector(self.elSet,
                                                                 gaussPts, 
                                                                 definition['result'])
            
#            if 'gaussPt' in definition and not isInteger(definition['gaussPt'] ):
#                # results are requested for multiple gausspoints
#                # this means that for every gaussPt a request is handed to each element,
#                # and the resulting array has one additional dimension ([elements, gaussPts, results])
#                gpts  = strToRange(definition['gaussPt'])
#                self.permanentElResultMemory = [ [el.getPermanentResultPtr(** dict(requestDictForElement, **{'gaussPt':gpt}) ) for gpt in gpts ] for el in self.elSet]
#            
#            else:
#                # either a single gaussPt, or no gaussPt at all was requested.
#                # Anyways, only a single request has to be handed to the elements
#                if 'gaussPt' in definition:
#                    requestDictForElement['gaussPt'] = int(definition['gaussPt'])
#                self.permanentElResultMemory = [el.getPermanentResultPtr(**requestDictForElement) for el in self.elSet]

        if 'f(x)' in definition:
            self.f = createMathExpression( definition['f(x)'] )
        else:
            self.f = None
            
        self.timeHistory = []
        
        # handle export of the fieldout at the end of the job:
        self.export = definition.get('export', False)
        if 'f_export(x)' in definition:
            self.f_export =  createMathExpression( definition['f_export(x)'] ) 
        else:
            self.f_export = None
        
    def getLastResult(self, **kw):
        return self.result[-1] if self.appendResults else self.result
    
    def getResultHistory(self, ):
        if not self.appendResults:
            raise Exception('fieldOuput {:} does not save any history; please define it with saveHistory=True!'.format(self.name))
        return np.asarray(self.result) 
    
    def getTimeHistory(self, ):
        return np.asarray( self.timeHistory )
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        
        tic = getCurrentTime()
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        self.timeHistory.append ( totalTime + dT )
        
        incrementResult = None
        
        if self.type == 'perNode':
            resVec = U if self.resultVector == 'U' else P
            incrementResult =  resVec [ self.node.fields[ self.field ] ]
            
        elif self.type == 'perNodeSet':
            resVec = U if self.resultVector == 'U' else P
#            incrementResult =  np.array( [  resVec [ n.fields[ self.field ] ] for n in self.nSet])
            incrementResult =  resVec[self.resultIndices].reshape( len(self.nSet), -1)
        
        elif self.type == 'perElementSet' or self.type == 'perElement':
#            incrementResult = np.asarray( self.permanentElResultMemory )
            incrementResult = self.elementResultCollector.getCurrentResults()
            
        if self.f:
            incrementResult = self.f( incrementResult )
        
        if self.appendResults:
            self.result.append(incrementResult)
        else:
            self.result = incrementResult
        
        toc = getCurrentTime()
        self.timeTotal += toc - tic
        
    def finalizeStep(self, U, P,):
        pass
    
    def finalizeJob(self, U, P,):
        
        print(self.name + "   " + str(self.timeTotal))
        
        if self.export:
            res = np.asarray ( self.result )
            if self.f_export:
                res = self.f_export(res)
            if res.ndim > 2:
                self.journal.message('Reshaping fieldOutput result for export in .csv file', self.name)
                res = res .reshape ( ( res.shape[0] , -1) )
            if self.appendResults and res.shape[0] == len(self.timeHistory): 
                # we also store the time, if result shape and time history are 'compatible'
                self.journal.message('Adding time history for export in .csv file', self.name)
                time = np.asarray ( self.timeHistory ).reshape(-1, 1)
                resultTable = np.hstack ( [time , res]   )
            else:
                resultTable = res
            
            np.savetxt('{:}.csv'.format( self.export ), resultTable, )
            

class FieldOutputController:
    """
    The central module for managing field outputs, which can be used by output managers
    """
    
    def __init__(self, modelInfo, inputFile, journal):
        
        self.fieldOutputs = {}
        
        if not inputFile['*fieldOutput']:
            return
        definition = inputFile['*fieldOutput'][0]
        
        for defLine in definition['data']:
            fpDef = stringDict( defLine ) 
            self.fieldOutputs [ fpDef['name'] ] = FieldOutput ( modelInfo, fpDef , journal)
        
    def finalizeIncrement(self, U, P, increment):
        for output in self.fieldOutputs.values():
            output.finalizeIncrement(U,P,increment)
            
    def finalizeStep(self, U, P,):
        for output in self.fieldOutputs.values():
            output.finalizeStep(U,P)
            
    def initializeStep(self, step, stepActions, stepOptions):
        for output in self.fieldOutputs.values():
            output.initializeStep( step, stepActions, stepOptions)
    
    
    def finalizeJob(self, U, P,):
        for output in self.fieldOutputs.values():
            output.finalizeJob(U,P)
        
    def getRequestData(self, request):
        pass
