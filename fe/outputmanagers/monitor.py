#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 21:26:01 2017

@author: matthias
"""

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
            
    def finalizeStep(self, U, P):
        pass
    
    def finalizeJob(self,U, P):
        pass