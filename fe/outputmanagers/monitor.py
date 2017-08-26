#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 21:26:01 2017

@author: matthias

A simple monitor to check results in the console during analysis.

Datalines:
"""
documentation = {'fieldOutput' : 'fieldOutput to be monitored',
                 'f(x)' : '(optional), apply math on the increment fieldOutput',
                 }

from fe.outputmanagers.outputmanagerbase import OutputManagerBase
from fe.utils.misc import stringDict
import sympy as sp

class OutputManager(OutputManagerBase):
    """ Simple monitor for nodes, nodeSets, elements and elementSets """
    
    identification = "Monitor"
    printTemplate = "{:}, {:}: {:}"
    
    def __init__(self, name, definitionLines, jobInfo, modelInfo, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []
        self.fieldOutputController = fieldOutputController
        
        for defline in definitionLines:
            entry = {}
            defDict = stringDict(defline)
            entry['fieldOutput'] = fieldOutputController.fieldOutputs [ defDict['fieldOutput'] ]
            f = defDict.get('f(x)', 'x')
            entry['f(x)'] = sp.lambdify ( sp.DeferredVector('x'), f , 'numpy')
#            entry['export'] = defDict.get('export', False)
#            if entry['export']:
#                entry['history'] = []
#            
            self.monitorJobs.append(entry)
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        for nJob in self.monitorJobs:
            result = nJob['f(x)'] ( nJob['fieldOutput'].getLastResult() )
            self.journal.message(self.printTemplate.format(nJob['fieldOutput'].name, 
                                                           nJob['fieldOutput'].type,
                                                           result),
                                 self.identification)
    
#            if nJob['export']:
#                nJob['history'].append(result)
            
    def finalizeStep(self, U, P):
        pass
    
    def finalizeJob(self,U, P):
        pass
#        exportfiles = defaultdict(list)
#        
#        for nJob in self.monitorJobs:
#            if nJob['export']:
#                exportfiles[ nJob['export'] ] .append(nJob['history'])
#        
#        for exportName, exportTable in exportfiles.items():
#            np.savetxt('{:}.csv'.format(exportName), np.asarray(exportTable).T)
