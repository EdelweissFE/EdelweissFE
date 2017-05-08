#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 11:37:26 2017

@author: matthias
"""

from fe.outputmanagers.outputmanagerbase import OutputManagerBase

from fe.utils.misc import stringDict
import numpy as np
import matplotlib.pyplot as plt

class OutputManager(OutputManagerBase):
    identification = "PathPlotter"
    
    def __init__(self, name, definitionLines, jobInfo, modelInfo, journal):
        self.journal = journal
        self.monitorJobs = []

        nodes = modelInfo['nodes']
        self.figures = {}
        
        for defline in definitionLines:
            entry = {}
            defDict = stringDict(defline)
            figName = defDict.get('figure', 'defaultFigure')
            
            if figName not in self.figures:
                self.figures[figName] = plt.figure()
            entry['figure'] = self.figures[figName]
            
            nSetName = entry['nSetName'] = defDict['nSet']
            nodes = modelInfo['nodeSets'][nSetName]
            field = entry['field'] = defDict['field']
            direct = entry['dir'] = int(defDict.get('direction', 1) )-1
            entry['result'] = defDict.get('result', 'U')
            entry['label'] = defDict.get('label', 'result along path')
            entry['resultIndices'] = [node.fields[field][direct] for node in nodes]
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
            location = U if nJob['result'] == 'U' else P
                
            indices = nJob['resultIndices']
            result =  location[indices]  
            
            if nJob['normalize']:
                result /= np.max(np.abs(result))
            
            figure = nJob ['figure'] 
            figure.gca().plot( nJob['pathDistances'], result, label=nJob['label'] )
        
        plt.legend()
        plt.show()
