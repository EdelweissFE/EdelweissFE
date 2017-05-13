#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:12:40 2017

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
        self.magnitudeAtStepStart = 0.0
        
        action = stringDict(definition)        
        self.surface = modelInfo['surfaces'][action['surface']]
        self.loadType = action['type']
        magnitude = np.asarray( [float(action['magnitude'])], dtype=np.double)
        
        self.delta = magnitude
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x
            
        self.idle = False
            
    def finishStep(self):
        self.magnitudeAtStepStart += self.delta
        self.idle = True
    
    def updateStepAction(self, definition, jobInfo, modelInfo, journal):
        action = stringDict(definition)
        self.delta = np.fromstring(action['magnitude'], dtype=np.double) - self.magnitudeAtStepStart 
        self.idle = False
    
    def getCurrentMagnitude(self, increment):
        if self.idle:
            t = 1.0
        else:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            t = stepProgress
        
        return self.magnitudeAtStepStart + self.delta * self.amplitude(t)
