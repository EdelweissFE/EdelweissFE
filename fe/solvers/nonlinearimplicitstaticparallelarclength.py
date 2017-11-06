#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:43:01 2017

@author: m9

(Parallel) Arc Length Solver, based on the proposed approach in Jirásek/Bažant 2001.
Replaces the NewtonRaphson scheme of the NISTParallel Solver.
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
        
        options =                   stepOptions['NISTArcLength']
        arcLengthControllerType =   options['arcLengthController']
        self.arcLengthController =  stepActions[ arcLengthControllerType ] [arcLengthControllerType] # name = module designation
        
        if 'stopCondition' in options:
             self.checkConditionalStop = createModelAccessibleFunction( options['stopCondition'], self.modelInfo, self.fieldOutputs  )
        else:
            self.checkConditionalStop = lambda : False
        
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
        """Arc length method to solve for an increment
        Implementation based on the proposed approach by """
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        
        if incNumber > 1 and self.checkConditionalStop():
            raise ConditionalStop
        
        iterationCounter =          0
        incrementResidualHistory =  dict.fromkeys( self.fieldIndices, (0.0, 0 ) )
        
        dirichlets =            activeStepActions['dirichlets']
        concentratedLoads =     activeStepActions['concentratedLoads']
        distributedDeadLoads =  activeStepActions['distributedDeadLoads']
        
        activeDeadLoads =  concentratedLoads or distributedDeadLoads
        
        R_ =            np.tile(P, (2,1)).T # 2 RHSs
        R_0 =           R_[:,0]
        R_f =           R_[:,1]
        F =             np.zeros_like(P)    # accumulated Flux vector 
        Pdeadloads =    np.zeros_like(P)
        Pdeadloads_f =  np.zeros_like(P)
        
        ddU = None
        
        Lambda =  self.Lambda
        dLambda = self.dLambda
        ddLambda = 0.0
        
        dU, isExtrapolatedIncrement, dLambda = self.extrapolateLastIncrement(extrapolation, increment, dU, dirichlets, lastIncrementSize, dLambda)
        
        referenceIncrement = incNumber, 1.0, 1.0, 0.0, 0.0, 0.0 
        
        Pdeadloads_f = self.assembleDeadLoads (Pdeadloads_f, concentratedLoads, distributedDeadLoads, I, referenceIncrement)
        
        while True:
            for geostatic in activeStepActions['geostatics']: geostatic.apply() 
            
            modifiedIncrement = incNumber, dLambda, Lambda + dLambda, 0.0, 0.0, 0.0

            P, V, F = self.computeElements(U, dU, P, V, I, J, F, increment)
            R_f[:] = 0.0
            
            if activeDeadLoads:
                Pdeadloads = self.assembleDeadLoads (Pdeadloads, concentratedLoads, distributedDeadLoads, I, modifiedIncrement)
                P   += Pdeadloads
                R_f += Pdeadloads_f
                
            R_0[:] = P
            
            R_0 = self.applyDirichlet( modifiedIncrement, R_0, dirichlets )
            R_f = self.applyDirichlet ( referenceIncrement, R_f, dirichlets)
            
            for constraint in self.constraints.values(): 
                R_0[constraint.globalDofIndices] = 0.0
                R_f[constraint.globalDofIndices] = 0.0
             
            if iterationCounter > 0:
                if self.checkConvergency(R_0, ddU, F, iterationCounter, incrementResidualHistory):
                    break
                
                if self.checkDivergingSolution (incrementResidualHistory , maxGrowingIter):
                    raise DivergingSolution('Residual grew {:} times, cutting back'.format(maxGrowingIter))
                
            if iterationCounter == maxIter:
                raise  ReachedMaxIterations("Reached max. iterations in current increment, cutting back")
            
            K = self.assembleStiffness(V, I, J)
            K = self.applyDirichletK(K, dirichlets)
            
            ddU_ = self.linearSolve(K, R_ )
            
            # solve 2 eq. systems at once:
            # q_0 = K⁻¹ * (  Pext_0  + dLambda * Pext_Ref - PInt  )
            # q_f = K⁻¹ * (  Pext_Ref  )
            ddU_0, ddU_f = ddU_[:,0], ddU_[:,1]
            
            # compute the increment of the load parameter. Method depends on the employed arc length controller
            ddLambda = self.arcLengthController.computeDDLambda( dU, ddU_0, ddU_f, increment ) 
            
            # assemble total solution
            ddU = ddU_0 + ddLambda * ddU_f
            
            dU +=       ddU
            dLambda +=  ddLambda
            
            iterationCounter += 1
           
        self.Lambda += dLambda
        self.dLambda = dLambda
        
        self.arcLengthController.finishIncrement(U, dU, dLambda)  
        
        return dU, iterationCounter, incrementResidualHistory
    
    def extrapolateLastIncrement(self, extrapolation, increment, dU, dirichlets, lastIncrementSize, dLambda):
        """ also account for dLambda
        seems to be not working properly yet"""
        
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        
        if extrapolation == 'linear' and lastIncrementSize:
            dLambda = dLambda * (incrementSize/lastIncrementSize) 
        else:
            dLambda = 0.0    
            
        dU, isExtrapolatedIncrement = super().extrapolateLastIncrement( extrapolation, increment, dU, dirichlets, lastIncrementSize)
        
        return dU, isExtrapolatedIncrement, dLambda