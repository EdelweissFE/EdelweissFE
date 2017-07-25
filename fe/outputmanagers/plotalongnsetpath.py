#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 11:37:26 2017

@author: matthias
"""

from fe.outputmanagers.outputmanagerbase import OutputManagerBase

from fe.utils.misc import stringDict
import numpy as np
import sympy as sp

class OutputManager(OutputManagerBase):
    identification = "PathPlotter"
    
    def __init__(self, name, definitionLines, jobInfo, modelInfo, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []

        nodes = modelInfo['nodes']
        self.plotter = plotter
        
        for defline in definitionLines:
            entry = {}
            defDict = stringDict(defline)
            entry['fieldOutput'] = fieldOutputController.fieldOutputs [ defDict['fieldOutput'] ]
            nodes = entry['fieldOutput'].nSet
            
            f = defDict.get('f(x)', 'x')
            entry['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), f , 'numpy')
            entry['label'] = defDict.get('label', entry['fieldOutput'].name)

            entry['figure'] = defDict.get('figure', 1)
            entry['axSpec'] = defDict.get('axSpec', 111)
            
            entry['normalize'] = defDict.get('normalize', False)
            self.monitorJobs.append(entry)
            
            #compute distance(s), node 0 is the reference node in the 'origin'
            entry['pathDistances'] = [0.0]
            # 1) distances between nodes:
            distancesBetweenNodes = [ np.linalg.norm( nodes[i+1].coordinates - nodes[i].coordinates )  for i in range( len(nodes) - 1 ) ]
            # 2) distances with respect to node 0
            for dist in distancesBetweenNodes :
                entry['pathDistances'].append( entry['pathDistances'][-1] + dist)
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        pass
            
    def finalizeStep(self, U, P,):
        pass
    
    def finalizeJob(self, U, P,):
        for nJob in self.monitorJobs:

            result = nJob['fieldOutput'].getLastResult()
            result = nJob['f(x)'] ( result )
                    
            if result.ndim >= 2:
                raise Exception('plot along nset path: result ndim >=2 ')

            if nJob['normalize']:
                result /= np.max(np.abs(result))
            
            self.plotter.plotXYData(nJob['pathDistances'], result, 
                                    nJob['figure'], nJob['axSpec'], nJob )
