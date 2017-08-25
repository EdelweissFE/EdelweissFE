#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 14:57:48 2017

@author: matthias
"""

import numpy as np
import sympy as sp

from fe.utils.misc import stringDict
from fe.utils.meshtools import  extractNodesFromElementSet

class FieldOutput:
    """
    Entity of a fieldOutput request
    """
    def __init__(self, modelInfo, definition, journal):
        self.name = definition['name']
        self.journal = journal
        
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
            self.element = modelInfo['elements'] [ int(definition['element']) ]
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
                
        elif self.type == 'perElement' or self.type == 'perElementSet':
            
            requestDictForElement = {}
            requestDictForElement['result'] = definition['result']
            if 'index' in definition:
                idcs = definition['index']
                if ':' in definition['index']:
                    idcs=[int (i) for i in idcs.split(':')]
                    requestDictForElement['idxStart'], requestDictForElement['idxStop'] = idcs
                else:
                    idx = int (idcs)
                    requestDictForElement['idxStart'], requestDictForElement['idxStop'] = idx, idx+1
            if 'gaussPt' in definition:
                requestDictForElement['gaussPt'] = int( definition['gaussPt'] )
            
            if self.type == 'perElement':
                self.permanentElResultMemory = [ self.element.getPermanentResultPtr(**requestDictForElement) ]
            elif self.type == 'perElementSet':
                self.permanentElResultMemory = [el.getPermanentResultPtr(**requestDictForElement) for el in self.elSet]
            
        if 'f(x)' in definition:
            self.fString = definition['f(x)']
            self.f = sp.lambdify ( sp.DeferredVector('x'), definition['f(x)'] , 'numpy')
        else:
            self.f = None
            
        self.timeHistory = []
        
    def getLastResult(self, **kw):
        return self.result[-1] if self.appendResults else self.result
    
    def getResultHistory(self, ):
        if not self.appendResults:
            raise Exception('fieldOuput {:} does not save any history; please define saveHistory=True!'.format(self.name))
        return np.asarray(self.result) 
    
    def getTimeHistory(self, ):
        return np.asarray( self.timeHistory )
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        self.timeHistory.append ( totalTime + dT )
        
        incrementResult = None
        
        if self.type == 'perNode':
            resVec = U if self.resultVector == 'U' else P
            incrementResult =  resVec [ self.node.fields[ self.field ] ]
            
        elif self.type == 'perNodeSet':
            resVec = U if self.resultVector == 'U' else P
            incrementResult =  np.array( [  resVec [ n.fields[ self.field ] ] for n in self.nSet])
        
        elif self.type == 'perElementSet' or self.type == 'perElement':
            incrementResult = np.asarray( self.permanentElResultMemory )
            
        if self.f:
            incrementResult = self.f( incrementResult )
        
        if self.appendResults:
            self.result.append(incrementResult)
        else:
            self.result = incrementResult
            
    def finalizeStep(self, U, P,):
        pass
    
    def finalizeJob(self, U, P,):
        pass
    

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
