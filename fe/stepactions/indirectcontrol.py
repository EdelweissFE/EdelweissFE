#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:35:44 2017

@author: matthias
"""

documentation={
        
        'nSet':'nSet for application of the BC',
        '1,2,3':'prescribed values in directions',
        'field': 'field for BC',
        'f(t)':'(optional) define an amplitude',
        
        }

from fe.stepactions.stepactionbase import StepActionBase
import numpy as np

class StepAction(StepActionBase):
    """ Dirichlet boundary condition, based on a node set """
    def __init__(self, name, action, jobInfo, modelInfo, journal):
                
        self.name = name
        self.journal = journal
        
        if 'deactive' in action:
            self.active = False
            return
        else:
            self.active = True
        
        nodes = modelInfo['nodes']
        nodeSets = modelInfo['nodeSets']
        
        controltype = action['controltype']
        
        if controltype == "crackmouthopening":
            self.c = np.array([-1, 1])
        
        self.L =  float ( action['L'] )
        
        self.dof1 = eval ( action['dof1'] )
        self.dof2 = eval ( action['dof2'] )
        
        
        self.idcs = np.array([self.dof1, self.dof2])
        
    def computeDDLambda(self, dU, ddU_0, ddU_f, increment ):
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        dL = incrementSize * self.L
        
        ddLambda = ( dL - self.c.dot ( dU [self.idcs]  + ddU_0 [self.idcs] ) )  /  self.c.dot ( ddU_f [self.idcs] )
        return ddLambda
    
    def finishIncrement(self, U, dU, dLambda):
        self.journal.message('Dof 1: {:5.5f}, Dof 2: {:5.5f}'.format( 
                U [self.dof1] + dU [self.dof1],
                U [self.dof2] + dU [self.dof2]), self.name )
        
    def finishStep(self,):
        pass
#        self.active = False
    
    def updateStepAction(self, action):
        pass