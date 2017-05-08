#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 13:03:09 2017

@author: matthias
"""
from fe.stepactions.stepactionbase import StepActionBase
from fe.utils.misc import stringDict
import numpy as np
import sympy as sp

class StepAction(StepActionBase):
    def __init__(self, name, definition, jobInfo, modelInfo, journal):
        """ create dirichlet dictionary with node boundary condition in 
        keytype 'indices': array of global dof indices
                'delta':   prescribed deltaValue """
                
        self.name = name
        dirichletIndices = []
        dirichletDelta = []
        
        nodeSets = modelInfo['nodeSets']
        action = stringDict(definition) 
        field = action['field']
        
        for x, direction  in enumerate(['1', '2', '3']):
            if direction in action:
                directionIndices = [node.fields[field][x] for node in nodeSets[action['nSet']]]
                dirichletIndices += directionIndices
                dirichletDelta += [float(action[direction])] * len(directionIndices)
                            
        self.indices = np.array(dirichletIndices)
        self.delta = np.array(dirichletDelta)
        
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x
        
    
    def updateStepAction(self, definitionLines, jobInfo, modelInfo, journal):
        pass
    
    def getDelta(self, increment):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        return self.delta * ( self.amplitude ( stepProgress ) - 
                             (self.amplitude ( stepProgress - incrementSize )))
