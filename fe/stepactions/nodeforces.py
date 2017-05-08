#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 19:33:06 2017

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
        nodeForceIndices = []
        nodeForceDelta = []
        
        nodeSets = modelInfo['nodeSets']
        
        action = stringDict(definition)        
        field = action['field']
        for x, direction  in enumerate(['1', '2', '3']):
            if direction in action:
                directionIndices = [node.fields[field][x] for node in nodeSets[action['nSet']]]
                nodeForceIndices += directionIndices
                nodeForceDelta += [float(action[direction])] * len(directionIndices)
                            
        self.indices = np.array(nodeForceIndices)
        self.deltaP = np.array(nodeForceDelta)
        
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x
        
    
    def updateStepAction(self, definitionLines, jobInfo, modelInfo, journal):
        pass
    
    def applyOnP(self, P, increment):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        P[self.indices] += self.deltaP * self.amplitude( stepProgress )
        return P