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

class OutputManager(OutputManagerBase):
    identification = "NodeSetMonitor"
    printTemplate = "nSet {:}, {:} {:}: {:}"
    
    def __init__(self, name, definitionLines, jobInfo, modelInfo, journal):
        self.journal = journal
        self.monitorJobs = []

        nodes = modelInfo['nodes']
        
        for defline in definitionLines:
            entry = {}
            defDict = stringDict(defline)
            nSetName = entry['nSetName'] = defDict['nSet']
            nodes = modelInfo['nodeSets'][nSetName]
            field = entry['field'] = defDict['field']
            direct = entry['dir'] = int(defDict['direction'])-1
            entry['result'] = defDict.get('result', 'U')
            entry['resultIndices'] = [node.fields[field][direct] for node in nodes]
            entry['export'] = defDict.get('export', False)
            
            f = defDict.get('f(x)', 'x')
            entry['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), f , 'numpy')
            
            if entry['export']:
                entry['history'] = []
            self.monitorJobs.append(entry)
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        for nJob in self.monitorJobs:
            
            location = U if nJob['result'] == 'U' else P
                
            indices = nJob['resultIndices']
            result = nJob['f(x)'] ( location[indices]  )
            self.journal.message(self.printTemplate.format(nJob['nSetName'], 
                                                           nJob['field'],  
                                                           nJob['result'], 
                                                           result),
                                 self.identification)
            if nJob['export']:
                nJob['history'].append(result)    
            
    def finalizeStep(self, U, P,):
        pass
    
    def finalizeJob(self, U, P,):
        exportfiles = defaultdict(list)
        
        for nJob in self.monitorJobs:
            if nJob['export']:
                exportfiles[ nJob['export'] ] .append(nJob['history'])
                
        for exportName, exportTable in exportfiles.items():
            np.savetxt('{:}.csv'.format(exportName), np.asarray(exportTable).T)