#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:35:44 2017

@author: matthias

Indirect (displacement) controller for the NISTArcLength solver
"""

documentation={
        
        'dof1':'Degree of freedom for the constraint ( model access expression )',
        'dof2':'Degree of freedom for the constraint ( model access expression )',
        'L' : 'Final distance (e.g. crack opening)',
        }

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np
from fe.utils.math import evalModelAccessibleExpression

class StepAction(StepActionBase):
    identification = 'IndirectControl'
    
    def __init__(self, name, action, jobInfo, modelInfo, journal):
                
        self.name =     name
        self.journal = journal
        self.modelInfo = modelInfo 
        self.c = np.array([-1, 1])
        self.currentL0 = 0.0
        
        self.L =  float ( action['L'] )
        self.dof1 = evalModelAccessibleExpression ( action['dof1'], modelInfo, fieldOutputs=None)
        self.dof2 = evalModelAccessibleExpression ( action['dof2'], modelInfo, fieldOutputs=None)
       
        self.definition = str( action.get('definition', 'absolute') ) 
        self.idcs = np.array([self.dof1, self.dof2])
        
    def computeDDLambda(self, dU, ddU_0, ddU_f, increment ):
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        dL = incrementSize * self.L
       
        ddLambda = ( dL - self.c.dot ( dU [self.idcs]  + ddU_0 [self.idcs] ) )  /  self.c.dot ( ddU_f [self.idcs] )
        return ddLambda
    
    def finishIncrement(self, U, dU, dLambda):
        self.journal.message('Dof 1: {:5.5f}, Dof 2: {:5.5f}'.format( 
                U [self.dof1] + dU [self.dof1],
                U [self.dof2] + dU [self.dof2]), self.identification )
        
    def finishStep(self, U, P):
        self.currentL0 =  self.c.dot ( U [self.idcs]  )
    
    def updateStepAction(self, action):
        if self.definition=='absolute':
            self.L =  float ( action['L'] ) - self.currentL0
        else:
            self.L =  float ( action['L'] )
        self.dof1 = evalModelAccessibleExpression ( action['dof1'], self.modelInfo, fieldOutputs=None)
        self.dof2 = evalModelAccessibleExpression ( action['dof2'], self.modelInfo, fieldOutputs=None)
        
        self.idcs = np.array([self.dof1, self.dof2])
        self.c = np.array([-1, 1])
