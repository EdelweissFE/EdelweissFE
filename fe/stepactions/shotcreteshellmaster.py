#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 14:10:20 2017

@author: matthias

Module for applying displacements on a shotcrete shell
"""

documentation={
        
        'nSet':'nSet for application of the BC',
        }

from fe.stepactions.stepactionbase import StepActionBase
from fe.utils.misc import stringDict
import numpy as np

class StepAction(StepActionBase):
    """ Dirichlet boundary condition, based on a node set """
    def __init__(self, name, definition, jobInfo, modelInfo, journal):
                
        self.name = name
        
        dirichletIndices = []
        
        nodeSets = modelInfo['nodeSets']
        action = stringDict(definition) 
        
        displacementsFile = action['displacements']
        nSet = action['nSet']
        
        x = np.loadtxt( displacementsFile )
        
        self.t = x[:,0]
        self.t[:] -= self.t[0] 
        
        self.U = x[:,1:]

        nodes =  nodeSets [ nSet ]
        
        dirichletIndices = [node.fields['displacement'] for node in nodes]
        self.indices = np.array(dirichletIndices).ravel()
  
        
    def finishStep(self,):
        pass
    
    def updateStepAction(self, definition):
        pass

    def getDelta(self, increment):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        if dT == 0.0:
            return np.zeros_like(self.indices)
        
        delta = np.array([np.interp( (totalTime + dT),  self.t  , x ) - 
                          np.interp( (totalTime     ),  self.t  , x ) for x in self.U.T] )
        
        return delta