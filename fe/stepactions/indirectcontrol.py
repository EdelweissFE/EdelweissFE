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
        
#        n1 = nodes[ int ( action['node1'] )  ]
#        n2 = nodes[ int ( action['node2'] ) ]
#        
#        field = action['field']
        self.L =  float ( action['L'] )
        
#        dof1 = n1.fields[field][ int (action['dir1']) -1 ]
#        dof2 = n2.fields[field][ int (action['dir2']) -1 ]
        
        dof1 = eval ( action['dof1'] )
        dof2 = eval ( action['dof2'] )
        
#        if 'stopAt' in action:
            
#            self.stop
            
        
        self.idcs = np.array([dof1, dof2])
        
    def computeDDLambda(self, dU, ddU_0, ddU_f, increment ):
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        dL = incrementSize * self.L
        
        ddLambda = ( dL - self.c.dot ( dU [self.idcs]  + ddU_0 [self.idcs] ) )  /  self.c.dot ( ddU_f [self.idcs] )
        return ddLambda
    
    def finishIncrement(self,):
#        self.journal.message('node 1: {:}, node 2: {:}, difference: {:}'.format( , self.name )
        pass
        
    def finishStep(self,):
        pass
#        self.active = False
    
    def updateStepAction(self, action):
        pass