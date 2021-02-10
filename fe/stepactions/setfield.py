#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 9 10:05:41 2021

@author: matthias
"""

documentation={
        'fieldOutput':'field output to set',
        'type':'const',
        'value':'values',
        }

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np

class StepAction(StepActionBase):
    """ Defines node based load, defined on a nodeset."""
    def __init__(self, name, action, jobInfo, modelInfo, fieldOutputController, journal):

        self.name = name
        self.active = True
        self.journal = journal
        self.fieldOutputName = action['fieldOutput']
        self.fieldOutput = fieldOutputController.fieldOutputs [ self.fieldOutputName ]
        self.type = action['type']
        self.value = action['value']
        
    def finishStep(self, U, P, stepMagnitude=None):
        self.active= False

    def updateStepAction(self, action):
        self.active= True

    def apply(self, ):
        if self.type == 'uniform':
            currentResults = np.zeros_like( self.fieldOutput.getLastResult() )
            currentResults[:] = float ( self.value )
            self.journal.message('setting field {:} to uniform value {:}'.format(self.fieldOutputName, self.value),  self.identification )
            self.fieldOutput.setResults ( currentResults )
        
