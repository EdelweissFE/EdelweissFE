#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 21:10:50 2017

@author: matthias
"""
from fe.outputmanagers.outputmanagerbase import OutputManagerBase
from fe.utils.misc import stringDict
from collections import defaultdict
import numpy as np
import sympy as sp

class OutputManager(OutputManagerBase):
    
    """ Simple monitor for nodes"""
    
    identification = "NodeMonitor"
    
    printTemplate = "node {:}, {:} {:}: {:}"
    
    def __init__(self, name, definitionLines, jobInfo, modelInfo, journal):
        self.journal = journal
        self.monitorJobs = []
        
        nodes = modelInfo['nodes']
        
        for defline in definitionLines:
            entry = {}
            defDict = stringDict(defline)
            entry['nodeNum'] = int(defDict['node'])
            entry['field'] = defDict['field']
            entry['result'] = defDict.get('result', 'U')
            entry['resultIndices'] = nodes[entry['nodeNum']].fields[defDict['field']]
            entry['export'] = defDict.get('export', False)
            if entry['export']:
                entry['history'] = []
                
            f = defDict.get('f(x)', 'x')
            entry['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), f , 'numpy')
            
            self.monitorJobs.append(entry)
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        for nJob in self.monitorJobs:
            
            location = U if nJob['result'] == 'U' else P
                
            indices = nJob['resultIndices']
            result = nJob['f(x)'] (location[indices])
            self.journal.message(self.printTemplate.format(nJob['nodeNum'], 
                                                           nJob['field'], 
                                                           nJob['result'],
                                                           result),
                                 self.identification)
            if nJob['export']:
                nJob['history'].append(result)    
            
    def finalizeStep(self, U, P):
        pass
    
    def finalizeJob(self,U, P):
        exportfiles = defaultdict(list)
        
        for nJob in self.monitorJobs:
            if nJob['export']:
                exportfiles[ nJob['export'] ] .append(nJob['history'])
                
        for exportName, exportTable in exportfiles.items():
            np.savetxt('{:}.csv'.format(exportName), np.asarray(exportTable).T)
    
    