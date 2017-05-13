#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 19:52:53 2017

@author: matthias
"""
from fe.stepactions.stepactionbase import StepActionBase
from fe.utils.misc import stringDict
import numpy as np


class StepAction(StepActionBase):
    def __init__(self, name, definition, jobInfo, modelInfo, journal):
        """ create dirichlet dictionary with node boundary condition in 
        keytype 'indices': array of global dof indices
                'delta':   prescribed deltaValue """
                
        self.name = name
        action = stringDict(definition) 
        
        self.geostaticElements = modelInfo['elSets'] [ action.get('elSet', 'all')]
        self.p1 = float(action['p1'])
        self.p2 = float(action.get('p2', 0.0))
        self.level1 = float(action.get('h1', 1.0))
        self.level2 = float(action.get('h2', -1.0))
        self.xLateral = float(action.get('xLateral', 1.0))
        self.yLateral = float(action.get('yLateral', 1.0))
        
        self.geostaticDefinition = np.array([
                self.p1,
                self.p2,
                self.level1,
                self.level2,
                self.xLateral,
                self.yLateral,
                ])
    
        self.active = True    
    
    def finishStep(self):
        self.active = False
    
    def updateStepAction(self, definitionLines, jobInfo, modelInfo, journal):
        pass
    
    
