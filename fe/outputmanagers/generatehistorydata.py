#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:39:30 2017

@author: matthias
"""
from fe.outputmanagers.outputmanagerbase import OutputManagerBase

from fe.utils.misc import stringDict
from collections import defaultdict
import numpy as np
import sympy as sp
from itertools import chain

class OutputManager(OutputManagerBase):
    identification = "generateHistoryData"
    
    def __init__(self, name, definitionLines, jobInfo, modelInfo, fieldOutputController, journal, plotter):

        self.outputData = {}
        for defline in definitionLines:
#            print(defline)
            defDict = stringDict(defline)
#            print(defDict)
            entry = {}
            entry['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), defDict.get('f(x)', 'x') , 'numpy')
            
            if 'time' in defDict.keys():
                entry['history'] = []
            
            elif 'node' in defDict.keys():
                entry['node'] = int(defDict.get('node',1))
                entry['result'] = defDict.get('result','U')
                direct = entry['dir'] = int(defDict['direction'])-1
                field = defDict.get('field','displacement')
                entry['resultIndices'] = modelInfo['nodes'][entry['node']].fields[field][direct]     
                entry['history'] = []

            elif 'nSet' in defDict.keys():          
                nSetName = entry['nSet'] = defDict['nSet']
                nodes = modelInfo['nodeSets'][nSetName]
                field = defDict.get('field','displacement')
                entry['result'] = defDict.get('result', 'U')
                entry['resultIndices'] = [node.fields[field][ int(defDict['direction'])-1 ] for node in nodes]
                entry['history'] = []
                
            elif 'element' in defDict.keys():    
                entry['element']  = int(defDict.get('element','1'))
                entry['result'] = defDict.get('result', 'sdv')
                entry['index'] = int(defDict['index'])
                entry['idxStart'] = int(defDict['index'])
                entry['idxStop'] = int(defDict['index'])+1
                entry['gaussPt'] =  int(defDict.get('gaussPt', '1'))
                entry['object'] = modelInfo['elements'][entry['element']]
                entry['history'] = []

            elif 'elSet' in defDict.keys():   
                entry['elSet']  = defDict.get('elSet','all')
                entry['result'] = defDict.get('result', 'sdv')
                entry['index'] = int(defDict['index'])
                entry['idxStart'] = int(defDict['index'])
                entry['idxStop'] = int(defDict['index'])+1
                entry['gaussPt'] =  int(defDict.get('gaussPt', '1'))
                entry['object'] =  modelInfo['elementSets'][entry['elSet']]
                entry['history'] = []
            
            self.outputData[defDict['type']] = entry
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        
        for typeID, outputJob in self.outputData.items():
            if 'time' in outputJob.keys():
                result = increment[2] # write end of step time
                
            elif ('node' in outputJob.keys()) or ('nSet' in outputJob.keys()):
                location = U if outputJob['result'] == 'U' else P
                
                indices = outputJob['resultIndices']
                result = outputJob['f(x)'] ( location[indices]  )
            
            elif ('element' in outputJob.keys()):
                if outputJob['result'] == 'stress' or outputJob['result'] == 'strain':
                    result = outputJob['f(x)'] ( outputJob['object'].getPermanentResultPtr(**outputJob)[outputJob['index']] )
                else:
                    result = outputJob['f(x)'] ( outputJob['object'].getPermanentResultPtr(**outputJob))
 
            elif ('elSet' in outputJob.keys()):
                if outputJob['result'] == 'stress' or outputJob['result'] == 'strain':
                    result = outputJob['f(x)'] ( [obj.getPermanentResultPtr(**outputJob)[outputJob['index']] for obj in outputJob['object'] ] )
                else:
                    result = outputJob['f(x)'] (  [obj.getPermanentResultPtr(**outputJob)  for obj in outputJob['object'] ] )
                   
            outputJob['history'].append(result)    
            
    def finalizeStep(self, U, P,):
        pass
    
    def finalizeJob(self, U, P,):
        
        # generate jobOutput dictionary which can be accessed by paramIdentificiation.py
        try:
            self.jobOutput
        except AttributeError:
            self.jobOutput = defaultdict(list)
        
        for dataType, outputJob in self.outputData.items():
            self.jobOutput[dataType].append(outputJob['history'])
            self.jobOutput[dataType] = list(chain.from_iterable(self.jobOutput[dataType])) # flatten list