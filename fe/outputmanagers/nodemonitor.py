#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 21:10:50 2017

@author: matthias
"""
from fe.utils.misc import stringDict
from collections import defaultdict
import numpy as np

class OutputManager:
    identification = "NodeMonitor"
    printTemplate = "node {:}, {:} {:} {:}: {:}"
    resultFunctions = {'1' : lambda x: x[0],
                       '2' : lambda x: x[1],
                       '3' : lambda x: x[2],
                       'all' : lambda x: x,
                       'sum' : lambda x: np.sum(x),
                       'mean' : lambda x: np.mean(x),
                       'magnitude' : lambda x: np.linalg.norm(x),
                       }
    
    def __init__(self, definitionLines, jobInfo, modelInfo, journal):
        self.journal = journal
        self.monitorJobs = []
        
        nodes = modelInfo['nodes']
        
        for defline in definitionLines:
            entry = {}
            defDict = stringDict(defline)
            entry['nodeNum'] = int(defDict['node'])
            entry['field'] = defDict['field']
            entry['result'] = defDict['result']
            entry['type'] = defDict.get('type', 'U')
            entry['resultFun'] = self.resultFunctions[ defDict['result']] 
#            entry['resultIndices'] = nodes[entry['nodeNum']]['fields'][defDict['field']]
            entry['resultIndices'] = nodes[entry['nodeNum']].fields[defDict['field']]
            entry['export'] = defDict.get('export', False)
            if entry['export']:
                entry['history'] = []
            self.monitorJobs.append(entry)
    
    def finalizeIncrement(self, U, P, increment):
        for nJob in self.monitorJobs:
            if nJob['type'] == 'U':
                location = U
            else:
                location = P
                
            indices = nJob['resultIndices']
            result = nJob['resultFun'] (location[indices])
            self.journal.message(self.printTemplate.format(nJob['nodeNum'], 
                                                           nJob['field'], 
                                                           nJob['type'], 
                                                           nJob['result'],
                                                           result),
                                 self.identification)
            if nJob['export']:
                nJob['history'].append(result)    
            
    def finalizeStep(self,):
        pass
    
    def finalizeJob(self,):
        exportfiles = defaultdict(list)
        
        for nJob in self.monitorJobs:
            if nJob['export']:
                exportfiles[ nJob['export'] ] .append(nJob['history'])
                
        for exportName, exportTable in exportfiles.items():
            np.savetxt('{:}.csv'.format(exportName), np.asarray(exportTable).T)
    
    