#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:08:32 2017

@author: matthias
"""


from fe.outputmanagers.outputmanagerbase import OutputManagerBase

from fe.utils.misc import stringDict
import numpy as np

class OutputManager(OutputManagerBase):
    identification = "TimeMonitor"
    def __init__(self, name, definitionLines, jobInfo, modelInfo, journal):
        self.journal = journal
        self.monitorJobs = []
        defDict = stringDict(definitionLines[0])
        self.exportFile = defDict['export']
        self.timeVals = []
    
    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        self.timeVals.append(totalTime + dT)

    def finalizeStep(self, U, P,):
        pass
    
    def finalizeJob(self, U, P,):
        np.savetxt('{:}.csv'.format(self.exportFile), np.asarray( self.timeVals).T)