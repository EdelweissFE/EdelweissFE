#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 19:33:06 2017

@author: matthias

Apply simple node forces on a nSet.
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
    """ Defines node based load, defined on a nodeset."""
    def __init__(self, name, action, jobInfo, modelInfo, journal):
                
        self.name = name
        nodeForceIndices = []
        nodeForceDelta = []
        nodeSets = modelInfo['nodeSets']
        
        self.field = action['field']
        self.idle = False
        self.nSet = nodeSets[action['nSet']]
        
        for x, direction  in enumerate(['1', '2', '3']):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                nodeForceIndices += directionIndices
                nodeForceDelta += [float(action[direction])] * len(directionIndices)
                 
        self.indices = np.asarray(nodeForceIndices, dtype=np.int)
        self.nodeForcesStepStart = np.zeros_like(self.indices, dtype=np.double)
        self.nodeForcesDelta = np.asarray(nodeForceDelta)
        
        self.amplitude = self.getAmplitude(action)
    
    def getAmplitude(self, action):
        if 'f(t)' in action:
            t = sp.symbols('t')
            amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            amplitude = lambda x:x
            
        return amplitude
        
    def finishStep(self):
        self.idle = True
        self.nodeForcesStepStart += self.nodeForcesDelta * self.amplitude(1)
    
    def updateStepAction(self, definition):
        self.idle = False
        
        action = stringDict(definition)
        nodeForceDelta = []
        for x, direction  in enumerate(['1', '2', '3']):
            if direction in action:
                directionIndices = [node.fields[self.field][x] for node in self.nSet]
                nodeForceDelta += [float(action[direction])] * len(directionIndices)
                 
        self.nodeForcesDelta = np.asarray(nodeForceDelta)
        self.amplitude = self.getAmplitude(action)
        
    def applyOnP(self, P, increment):
        if self.idle:
            P[self.indices] += self.nodeForcesStepStart
        else:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            t = stepProgress
            amp = self.amplitude ( t )
            P[self.indices] += self.nodeForcesStepStart +  self.nodeForcesDelta * amp
            
        return P