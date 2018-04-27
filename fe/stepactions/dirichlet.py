#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 13:03:09 2017

@author: matthias

Dirichlet boundary condition.
If not modified in subsequent steps, the BC is held constant.
"""
documentation={
        
        'nSet':'nSet for application of the BC',
        '1,2,3':'prescribed values in directions',
        'field': 'field for BC',
        'f(t)':'(optional) define an amplitude',
        
        }

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np
import sympy as sp

class StepAction(StepActionBase):
    """ Dirichlet boundary condition, based on a node set """
    
    def __init__(self, name, action, jobInfo, modelInfo, journal):
                
        self.name = name
        
        dirichletIndices = []
        dirichletDelta = []
        
        nodeSets = modelInfo['nodeSets']
        self.field = action['field']

        self.action = action
        self.nSet = nodeSets [ action['nSet'] ]

        for x, direction  in enumerate(['1', '2', '3']):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                dirichletIndices += directionIndices
                dirichletDelta += [float(action[direction])] * len(directionIndices)
                            
        self.indices = np.array(dirichletIndices)
        self.delta = np.array(dirichletDelta)
        
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x
        
        self.active = True
        
    def finishStep(self, U, P):
        
        self.active = False
    
    def updateStepAction(self, action):
        
        self.active = True
        dirichletIndices = []
        dirichletDelta = []
        
        for x, direction  in enumerate(['1', '2', '3']):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                dirichletIndices += directionIndices
                dirichletDelta += [float(action[direction])] * len(directionIndices)
                            
        self.delta = np.array(dirichletDelta)
        
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x

    def getDelta(self, increment):
        
        if self.active:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            return self.delta * ( self.amplitude ( stepProgress ) - 
                                 (self.amplitude ( stepProgress - incrementSize )))
        else:
            return 0.0
