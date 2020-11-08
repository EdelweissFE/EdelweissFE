#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mo July 29 10:50:53 2019

@author: magdalena

"""

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np


documentation = {
        'sdv' : 'set number of state variable for initialization',
        'val' : 'set initial value of state variable',
        }

class StepAction(StepActionBase):
    """ Initializes elements of set with an Abaqus-like geostatic stress state.
    Is automatically deactivated at the end of the step."""
    def __init__(self, name, action, jobInfo, modelInfo, journal):

        self.name = name
        
        self.theElements= modelInfo['elementSets'] [ action.get('elSet', 'all')]
        self.active = True    
        self.emptyDef = np.array([0.0])
    
    def finishStep(self, U, P, stepMagnitude=None):
        self.active = False
    
    def updateStepAction(self, definition):
        pass 
    
    def apply(self, ):
        for el in self.theElements:
            el.setInitialCondition('initialize material', self.emptyDef)
    
    
