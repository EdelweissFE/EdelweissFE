#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 13:15:14 2018

@author: matthias

Body force load.
If not modified in subsequent steps, the load held constant.
"""

documentation={
        'forceVector': 'force vector',
        'delta': 'in subsequent steps only: define the new force vector incrementally',
        'f(t)':'(optional) define an amplitude',
        }

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np
import sympy as sp

class StepAction(StepActionBase):
    """ Simple Body Force """
    def __init__(self, name, action, jobInfo, modelInfo, journal):
                
        self.name = name
        self.forceAtStepStart = 0.0
        self.elements = modelInfo['elementSets'] [ action['elSet'] ]
        magnitude = np.fromstring(action['forceVector'], sep=',', dtype=np.double)
        
        if len(magnitude) < jobInfo['domainSize']:
            raise Exception('BodyForce {:}: force vector has wrong dimension!'.format(self.name))
        
        self.delta = magnitude
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x
            
        self.idle = False
            
    def finishStep(self, U, P, stepMagnitude=None):
        
        if not self.idle:
            if stepMagnitude == None:
                # standard case
                self.forceAtStepStart += self.delta * self.amplitude(1.0)
            else:
                # set the 'actual' increment manually, e.g. for arc length method
                self.forceAtStepStart += self.delta * stepMagnitude
                
            self.delta = 0
            self.idle = True
    
    def updateStepAction(self, action):
        
        if 'forceVector' in action:
            self.delta =  np.fromstring(action['forceVector'], sep=',', dtype=np.double) - self.forceAtStepStart 
        elif 'delta' in action:
            self.delta = np.fromstring(action['delta'], sep=',', dtype=np.double)
            
        if 'f(t)' in action:
            t = sp.symbols('t')
            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
        else:
            self.amplitude = lambda x:x
            
        self.idle = False
    
    def getCurrentBodyForce(self, increment):
        if self.idle == True:
            t = 1.0
        else:
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            t = stepProgress
            
        return self.forceAtStepStart + self.delta * self.amplitude(t)