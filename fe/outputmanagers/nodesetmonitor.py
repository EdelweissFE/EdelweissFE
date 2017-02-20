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

class OutputManager(OutputManagerBase):
    identification = "NodeSetMonitor"
    printTemplate = "nSet {:}, {:} {:}: {:}"
    resultFunctions = {
                       'sum' : lambda x: np.sum(x),
                       'mean' : lambda x: np.mean(x),
                       }
    
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
            entry['type'] = defDict.get('type', 'U')
            entry['resultFun'] = self.resultFunctions[ defDict['result']] 
            entry['resultIndices'] = [node.fields[field][direct] for node in nodes]
            entry['export'] = defDict.get('export', False)
            if entry['export']:
                entry['history'] = []
            self.monitorJobs.append(entry)
    
    def initializeStep(self, step, stepActions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        for nJob in self.monitorJobs:
            if nJob['type'] == 'U':
                location = U
            else:
                location = P
                
            indices = nJob['resultIndices']
            result = nJob['resultFun'] ( location[indices]  )
            self.journal.message(self.printTemplate.format(nJob['nSetName'], 
                                                           nJob['field'], 
                                                           nJob['type'], 
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