#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 14:10:20 2017

@author: matthias

Module for applying displacements on a shotcrete shell
"""

documentation={
        
        'nSet':'nSet for application of the BC',
#        '1,2,3':'prescribed values in directions',
#        'field': 'field for BC',
#        'f(t)':'(optional) define an amplitude',
        }

from fe.stepactions.stepactionbase import StepActionBase
from fe.utils.misc import stringDict
import numpy as np
#import sympy as sp
import os

class StepAction(StepActionBase):
    """ Dirichlet boundary condition, based on a node set """
    def __init__(self, name, definition, jobInfo, modelInfo, journal):
                
        self.name = name
        
        dirichletIndices = []
#        dirichletDelta = []
        
        nodeSets = modelInfo['nodeSets']
        action = stringDict(definition) 
#        self.field = action['field']
        
#        displacementsFile = action['file']
        
        shotcreteMasterDirectory = action['directory']
        
        displacementsFile = np.loadtxt( os.path.join( shotcreteMasterDirectory, 'interpolatedDisplacements.csv' ) )
        
        x = displacementsFile
        
        self.t = x[:,0]
        self.t[:] -= self.t[0] 
        
        self.U = x[:,1:]

        self.action = action
        nodes =  nodeSets [ action['nSet'] +'Sorted']
        
#        for node in nodes:
            # assemble all indices of all nodes, on which the BCs should be imposed
        dirichletIndices = [node.fields['displacement'] for node in nodes]
        self.indices = np.array(dirichletIndices).ravel()
  
#        self.delta = np.zeros_like(dirichletDelta)
#        self.currentTotalDisplacement = np.zeros_like(dirichletDelta)
        
#        self.moving = True
        
    def finishStep(self,):
        pass
#        self.moving = False
    
    def updateStepAction(self, definition):
        pass
#        self.moving = True
#        action = stringDict(definition)
#
#        dirichletIndices = []
#        dirichletDelta = []
#        for x, direction  in enumerate(['1', '2', '3']):
#            if direction in action:
#                directionIndices = [node.fields[self.field][x] for node in self.nSet]
#                dirichletIndices += directionIndices
#                dirichletDelta += [float(action[direction])] * len(directionIndices)
#                            
#        self.delta = np.array(dirichletDelta)
#        
#        if 'f(t)' in action:
#            t = sp.symbols('t')
#            self.amplitude = sp.lambdify(t, sp.sympify(action['f(t)']), 'numpy')
#        else:
#            self.amplitude = lambda x:x

    def getDelta(self, increment):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        if dT == 0.0:
            return np.zeros_like(self.indices)
        
        delta = np.array([np.interp( (totalTime + dT),  self.t  , x ) - 
                          np.interp( (totalTime     ),  self.t  , x ) for x in self.U.T] )
        
        return delta
