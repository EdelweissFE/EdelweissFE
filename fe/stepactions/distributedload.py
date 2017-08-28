#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:12:40 2017

@author: matthias

Distributed load, applied on a surface set.
If not modified in subsequent steps, the load held constant.
"""

documentation={
        'surface':'surface for application of the distributed load',
        'magnitude':'dLoad magnitude',
        'delta': 'in subsequent steps only: define the new magnitude incrementally',
        'f(t)':'(optional) define an amplitude',
        }

from fe.stepactions.stepactionbase import StepActionBase
from fe.utils.misc import stringDict
import numpy as np
import sympy as sp

class StepAction(StepActionBase):
    """ Distributed load, defined on an element-based surface """
    def __init__(self, name, definition, jobInfo, modelInfo, journal):
                
        self.name = name
        self.magnitudeAtStepStart = 0.0
        
        action = stringDict(definition)        
        self.surface = modelInfo['surfaces'][action['surface']]
        self.loadType = action['type']
        magnitude = np.asarray([float(action['magnitude'])])
        
        self.delta = magnitude
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x
            
        self.idle = False
            
    def finishStep(self):
        self.magnitudeAtStepStart += self.delta * self.amplitude(1.0)
        self.delta=0
        self.idle = True
    
    def updateStepAction(self, definition):
        action = stringDict(definition)
        if 'magnitude' in action:
            self.delta = np.asarray([float(action['magnitude'])]) - self.magnitudeAtStepStart 
        elif 'delta' in action:
            self.delta = np.asarray([float(action['delta'])])   
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x
        self.idle = False
    
    def getCurrentMagnitude(self, increment):
        if self.idle == True:
            t = 1.0
        else:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            t = stepProgress
        return self.magnitudeAtStepStart + self.delta * self.amplitude(t)
