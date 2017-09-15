#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 19:52:53 2017

@author: matthias

@ Magdalena
Possible keywords for stepaction module geostatic
"""

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np


documentation = {
        'p1' : 'sig_x=sig_y=sig_z in first point',
        'h1' : 'y coordinate of first point, default=1.0',
        'p2' : 's11=s22=s33 in second point, default=p1',
        'h2' : 'y coordinate of second point, default=-1.0',
        'xLateral' : 'ratio of sig_x/sig_y, default=1.0',
        'zLateral' : 'ratio of sig_z/sig_y, default=1.0',
        }

class StepAction(StepActionBase):
    """ Initializes elements of set with an Abaqus-like geostatic stress state.
    Is automatically deactivated at the end of the step."""
    def __init__(self, name, action, jobInfo, modelInfo, journal):

        self.name = name
        
        self.geostaticElements = modelInfo['elementSets'] [ action.get('elSet', 'all')]
        self.p1 = float(action['p1'])
        self.p2 = float(action.get('p2', self.p1))
        self.level1 = float(action.get('h1', 1.0))
        self.level2 = float(action.get('h2', -1.0))
        self.xLateral = float(action.get('xLateral', 1.0))
        self.zLateral = float(action.get('zLateral', 1.0))
        
        self.geostaticDefinition = np.array([
                self.p1,
                self.level1,
                self.p2,
                self.level2,
                self.xLateral,
                self.zLateral,
                ])
    
        self.active = True    
    
    def finishStep(self):
        self.active = False
    
    def updateStepAction(self, definition):
        pass 
    
    def apply(self, ):
        for el in self.geostaticElements:
            el.setInitialCondition('geostatic stress', self.geostaticDefinition)
    
    
