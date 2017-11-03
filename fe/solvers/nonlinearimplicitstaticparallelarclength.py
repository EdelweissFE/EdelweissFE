#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:43:01 2017

@author: m9
"""

from fe.solvers.nonlinearimplicitstaticparallel import NISTParallel

import numpy as np
from fe.utils.exceptions import ReachedMaxIterations, DivergingSolution

class NISTPArcLength(NISTParallel):
    identification = "NISTPArcLength"
    
    
    def __init__(self, jobInfo, modelInfo, journal, fieldOutputController, outputmanagers):
        
        self.Lambda =  0.
        return super().__init__(jobInfo, modelInfo, journal, fieldOutputController, outputmanagers)
    
    def solveStep(self, step, time, stepActions, stepOptions, U, P):
        
        arcLengthControllerType = stepOptions['NISTArcLength']['arcLengthController']
        self.arcLengthController =  stepActions[ arcLengthControllerType ] [arcLengthControllerType] # name = module designation
        
        self.stepEndTime = 1.0
        self.totalEndTime = 1.0
        
        return super().solveStep(step, time, stepActions, stepOptions, U, P)
    
    def Newton(self, U, dU, 
           V, I, J, P, 
           dirichlets,
           concentratedLoads, 
           distributedDeadLoads, 
           activeGeostatics,
           increment,
           isExtrapolatedIncrement,
           maxIter,
           maxGrowingIter):
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        
        numberOfDofs =  self.nDof
        stepTimes = np.array([stepTime, totalTime])
        iterationCounter =          0
        incrementResidualHistory =  dict.fromkeys( self.fieldIndices, (0.0, 0 ) )
        activeDeadLoads =  concentratedLoads or distributedDeadLoads
        
        R_ =            np.tile(P, (2,1)).T
        R =             R_[:,0]
        R_f =           R_[:,1]
        F =             np.zeros_like(P)    # accumulated Flux vector 
        Pdeadloads =    np.zeros_like(P)
        
        ddU = None
        
        Lambda =  self.Lambda
        dLambda = 0.0
        ddLambda = 0.0
        
        referenceIncrement = incNumber, 1.0,  1.0,  1.0 * stepTime, self.stepEndTime, self.totalEndTime 
        
        while True:
            for geostatic in activeGeostatics: geostatic.apply() 
            
            P, V, F = self.computeElements(U, dU, stepTimes, dT, P, V, I, J, F)
            
            modifiedIncrement = incNumber, dLambda, Lambda + dLambda, 0.0, self.stepEndTime, self.totalEndTime
            
            if activeDeadLoads:
                Pdeadloads = self.assembleDeadLoads (Pdeadloads, concentratedLoads, distributedDeadLoads, I, stepTimes, dT, modifiedIncrement)
                P += Pdeadloads
                
            R[:] = P
            R_f[:] = 0.0
            
            R = self.applyDirichlet( modifiedIncrement, R, dirichlets)
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
            
            K = self.assembleStiffness(V, I, J, shape=(numberOfDofs, numberOfDofs) )
            K = self.applyDirichletK(K, dirichlets)
            
            ddU_ = self.linearSolve(K, R_ )
            
            ddU_0, ddU_f = ddU_[:,0], ddU_[:,1]
            
            ddLambda = self.arcLengthController.computeDDLambda( dU, ddU_0, ddU_f, increment ) 
            
            dU += ddU_0 + ddLambda * ddU_f
            dLambda += ddLambda
            
            iterationCounter += 1
           
#        self.arcLengthController.computeDDLambda( dU, ddU_0, ddU_f, increment )  
        self.Lambda += dLambda
        
        self.arcLengthController.finishIncrement()  
        return dU, iterationCounter, incrementResidualHistory