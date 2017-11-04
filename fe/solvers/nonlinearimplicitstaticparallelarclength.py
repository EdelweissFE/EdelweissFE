#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:43:01 2017

@author: m9
"""
from fe.solvers.nonlinearimplicitstaticparallel import NISTParallel

import numpy as np
from fe.utils.exceptions import ReachedMaxIterations, DivergingSolution, ConditionalStop
from fe.utils.math import createModelAccessibleFunction

class NISTPArcLength(NISTParallel):
    identification = "NISTPArcLength"
    
    def __init__(self, jobInfo, modelInfo, journal, fieldOutputController, outputmanagers):
        
        self.Lambda =  0.
        self.fieldOutputs = fieldOutputController.fieldOutputs
        self.modelInfo = modelInfo
        
        return super().__init__(jobInfo, modelInfo, journal, fieldOutputController, outputmanagers)
    
    def solveStep(self, step, time, stepActions, stepOptions, U, P):
        
        options = stepOptions['NISTArcLength']
        arcLengthControllerType = options['arcLengthController']
        
        if 'stopCondition' in options:
             self.checkConditionalStop = createModelAccessibleFunction( options['stopCondition'], self.modelInfo, self.fieldOutputs  )
        else:
            self.checkConditionalStop = False
        
        self.arcLengthController =  stepActions[ arcLengthControllerType ] [arcLengthControllerType] # name = module designation
        
        self.dLambda = 0.0
        
        return super().solveStep(step, time, stepActions, stepOptions, U, P)
    
    def solveIncrement(self, U, dU, 
           V, I, J, P, 
           activeStepActions,
           increment,
           lastIncrementSize,
           extrapolation,
           maxIter,
           maxGrowingIter):
        
        """Arc length method to solve for an increment"""
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        
        if self.checkConditionalStop and incNumber > 1.0 and self.checkConditionalStop():
            raise ConditionalStop
        
        iterationCounter =          0
        incrementResidualHistory =  dict.fromkeys( self.fieldIndices, (0.0, 0 ) )
        
        dirichlets =            activeStepActions['dirichlets']
        concentratedLoads =     activeStepActions['concentratedLoads']
        distributedDeadLoads =  activeStepActions['distributedDeadLoads']
        
        activeDeadLoads =  concentratedLoads or distributedDeadLoads
        
        R_ =            np.tile(P, (2,1)).T
        R =             R_[:,0]
        R_f =           R_[:,1]
        F =             np.zeros_like(P)    # accumulated Flux vector 
        Pdeadloads =    np.zeros_like(P)
        Pdeadloads_f =  np.zeros_like(P)
        
        ddU = None
        
        Lambda =  self.Lambda
        dLambda = self.dLambda
        ddLambda = 0.0
        
        dU, isExtrapolatedIncrement, dLambda = self.extrapolateLastIncrement(extrapolation, increment, dU, dirichlets, lastIncrementSize, dLambda)
        
        referenceIncrement = incNumber, 1.0,  1.0,  0.0, 0.0, 0.0 
        
        Pdeadloads_f = self.assembleDeadLoads (Pdeadloads_f, concentratedLoads, distributedDeadLoads, I, referenceIncrement)
        
        while True:
            for geostatic in activeStepActions['geostatics']: geostatic.apply() 
            
            modifiedIncrement = incNumber, dLambda, Lambda + dLambda, 0.0, 0.0, 0.0
            
            R_f[:] = 0.0
            P, V, F = self.computeElements(U, dU, P, V, I, J, F, increment)
            
            if activeDeadLoads:
                Pdeadloads = self.assembleDeadLoads (Pdeadloads, concentratedLoads, distributedDeadLoads, I, modifiedIncrement)
                P   += Pdeadloads
                R_f += Pdeadloads_f
                
            R[:] = P
            
            R = self.applyDirichlet( modifiedIncrement, R, dirichlets )
            R_f = self.applyDirichlet ( referenceIncrement, R_f, dirichlets)
            
            for constraint in self.constraints.values(): 
                R[constraint.globalDofIndices] = 0.0 # currently no external loads on rbs possible
             
            if iterationCounter > 0:
                if self.checkConvergency(R, ddU, F, iterationCounter, incrementResidualHistory):
                    break
                
                if self.checkDivergingSolution (incrementResidualHistory , maxGrowingIter):
                    raise DivergingSolution('Residual grew {:} times, cutting back'.format(maxGrowingIter))
                
            if iterationCounter == maxIter:
                raise  ReachedMaxIterations("Reached max. iterations in current increment, cutting back")
            
            K = self.assembleStiffness(V, I, J)
            K = self.applyDirichletK(K, dirichlets)
            
            ddU_ = self.linearSolve(K, R_ )
            
            ddU_0, ddU_f = ddU_[:,0], ddU_[:,1]
            
            ddLambda = self.arcLengthController.computeDDLambda( dU, ddU_0, ddU_f, increment ) 
            
            ddU = ddU_0 + ddLambda * ddU_f
            
            dU += ddU
            dLambda += ddLambda
            
            iterationCounter += 1
           
        self.Lambda += dLambda
        self.dLambda = dLambda
        
        self.arcLengthController.finishIncrement(U, dU, dLambda)  
        
        return dU, iterationCounter, incrementResidualHistory
    
    def extrapolateLastIncrement(self, extrapolation, increment, dU, dirichlets, lastIncrementSize, dLambda):
        """ also account for dLambda"""
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        
        if extrapolation == 'linear' and lastIncrementSize:
            dLambda = dLambda * (incrementSize/lastIncrementSize) 
        else:
            dLambda = 0.0    
            
        dU, isExtrapolatedIncrement = super().extrapolateLastIncrement( extrapolation, increment, dU, dirichlets, lastIncrementSize)
        
        return dU, isExtrapolatedIncrement, dLambda