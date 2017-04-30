#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 20:37:35 2017

@author: matthias
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from fe.utils.incrementgenerator import IncrementGenerator
from time import time as getCurrentTime

class NIST:
    """ This is the Nonlinear Implicit STatic -- solver.
    Designed to interface with Abaqus UELs
    Public methods are: __init__(), initializeUP() and solveStep(...).
    OutputManagers are updated at the end of each increment. """
    
    identification = "NISTSolver"
    
    defaultMaxInc = 1.0
    defaultMinInc = 1e-4
    defaultMaxNumInc = 1000
    defaultMaxIter = 10
    defaultCriticalIter = 5
    
    def __init__(self, jobInfo, modelInfo, journal, outputmanagers=None):
        self.nodes =        modelInfo['nodes']
        self.elements =     modelInfo['elements']
        self.nodeSets =     modelInfo['nodeSets']
        self.elementSets =  modelInfo['elementSets']
        self.fieldIndices = jobInfo['fieldIndices']
        
        self.flowCorrectionTolerances = jobInfo['flowCorrectionTolerance']
        self.effortResidualTolerances = jobInfo['effortResidualTolerance']
        self.effortResidualTolerancesAlt = jobInfo['effortResidualToleranceAlternative']
        
        # create headers for formatted output of solver
        nFields = len(self.fieldIndices.keys())
        self.iterationHeader = ("{:^25}"*nFields).format(*self.fieldIndices.keys())
        self.iterationHeader2 = (" {:<10}  {:<10}  ").format('||R||∞','||ddU||∞') *nFields
        self.iterationMessageTemplate = "{:11.2e}{:1}{:11.2e}{:1} "
        
        self.nDof = jobInfo['numberOfDofs']
        self.journal = journal
        self.outputmanagers = outputmanagers or []
        
        self.sizeVIJ = 0
        self.sizeNDofElementWise = 0
        
        for el in self.elements.values():
            self.sizeVIJ += el.sizeKe

        # create indices map to elements; V, I, J are of type np vectors
        # elementToIndexInVIJMap is a dictionary of {element : index in VIJ vectors} 
        V, I, J, elementToIndexInVIJMap = self.generateVIJ(self.elements, )
        self.V = V
        self.I = I
        self.J = J
        self.elementToIndexInVIJMap = elementToIndexInVIJMap # element  -> V[ .... idx ..  ]
        
    def initialize(self):
        """ Initialize the solver and return the 2 vectors for flow (U) and effort (P) """
        
        U = np.zeros(self.nDof)
        P = np.zeros(self.nDof)
        return U, P
        
    def solveStep(self, step, time, stepActions, U, P):
        """ Public interface to solve for an ABAQUS like step
        returns: boolean Success, U vector, P vector, and the new current total time """
            
        V = self.V
        J = self.J
        I = self.I
        numberOfDofs = self.nDof
        stepLength = step.get('stepLength', 1.0)
        extrapolation = stepActions['NISTSolverOptions'].get('extrapolation', 'linear')
        
        incGen = IncrementGenerator(time, 
                                    stepLength, 
                                    step.get('maxInc', self.defaultMaxInc), 
                                    step.get('minInc', self.defaultMinInc), 
                                    step.get('maxNumInc', self.defaultMaxNumInc), 
                                    self.journal)
        
        maxIter = step.get('maxIter', self.defaultMaxIter)
        criticalIter = step.get('crititcalIter', self.defaultCriticalIter)
        
        dU = np.zeros(numberOfDofs)
        Pext = np.zeros(numberOfDofs)
        
        # get indices where dirichlet BC are given
        dirichlet = stepActions['dirichlet']
        dirichletIndices = dirichlet['indices']
        
        # initialize pNewDt as np.array, for bypassing to external c++ code (UEL)
        pNewDT = np.zeros(1)
        lastIncrementSize = False
        
        computationTimeInElements = 0.0
        computationTimeInEqSystem = 0.0
        computationTimeMatrixConstruction = 0.0
        
        for increment in incGen.generateIncrement():
            incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
            
            self.journal.printSeperationLine()
            self.journal.message("increment {:}: {:8f}, {:8f}; time {:10f} to {:10f}".format(incNumber,
                                                                                        incrementSize, 
                                                                                        stepProgress,
                                                                                        totalTime,
                                                                                        totalTime + dT),
                                self.identification, level=1)
            self.journal.message(self.iterationHeader, self.identification, level=2)
            self.journal.message(self.iterationHeader2, self.identification, level=2)
            
            converged = False
            ddU = None
            iterationCounter = 0
            stepTimes = np.array([stepTime, totalTime])
            
            if extrapolation == 'linear' and lastIncrementSize:
                dU *= (incrementSize/lastIncrementSize) 
                dU = self.applyDirichlet(incrementSize, dU, dirichlet)      
                extrapolatedIncrement = True
            else:
                extrapolatedIncrement = False
                dU[:] = 0.0        

            while True:
                
                pNewDT[0] = 1e36
                
                tic =  getCurrentTime()
                P, V, pNewDT = self.computeElements(U, dU, stepTimes, dT, pNewDT, P, V, I, J,)
                toc =  getCurrentTime()
                computationTimeInElements += (toc-tic)
                
                if pNewDT[0] < 1.0:
                    self.journal.message("An element requests for a cutback", self.identification, level=2)
                    break
                    
                Pint  = P
                Pext[:] = 0.0
                
                if 'nodeForces' in stepActions:
                    # modify the load (effort) vector by external forces (effort)
                    nodeForces = stepActions['nodeForces']
                    Pext[nodeForces['indices']] += nodeForces['delta']
                    
                R = Pint - Pext
                
                
                
                if iterationCounter == 0 and not extrapolatedIncrement and dirichlet :
                    # first iteraion? apply dirichlet bcs and unconditionally solve
                    R = self.applyDirichlet(incrementSize, R, dirichlet)
                else:
                    # iteration cycle 1 or higher, time to check the convergency
                    R[dirichletIndices] = 0.0 # only entries not affected by dirichlet bcs contribute to the residual

                    converged = self.checkConvergency(R, ddU, iterationCounter)
                    
                if converged:
                    break
                
                if iterationCounter == maxIter-1:
                    self.journal.message("Reached max. iterations in current increment, cutting back", self.identification, 1)
                    break
                
                # not converged yet!
                # prepare global stiffness matrix
                
                tic =  getCurrentTime()
                K = coo_matrix( (V, (I,J)), shape=(numberOfDofs, numberOfDofs)).tocsr()
                toc =  getCurrentTime()
                computationTimeMatrixConstruction += toc - tic
                
                K_ = self.applyDirichletK(K, dirichlet)
                # solve !                
                
                tic =  getCurrentTime()
                ddU = spsolve(K_, R, )
                toc =  getCurrentTime()
                computationTimeInEqSystem += (toc-tic)
                dU += ddU
                
                iterationCounter += 1
                        
            if converged:
                U += dU
                lastIncrementSize = incrementSize
                if iterationCounter >= criticalIter:
                    incGen.preventIncrementIncrease()
                    
                for el in self.elements.values():
                    el.acceptLastState()
                self.journal.message("Converged in {:} iteration(s)".format(iterationCounter), self.identification, 1) 
                
                for man in self.outputmanagers:
                    man.finalizeIncrement(U, P, increment)
                
            else: 
                # get new increment by down-scaling of current increment
                if extrapolatedIncrement and iterationCounter == 0:
                    incGen.discardAndChangeIncrement(0.25)
                else:
#                    incGen.discardAndChangeIncrement(pNewDT[0] if pNewDT[0] < 1.0 else 0.25)
                    incGen.discardAndChangeIncrement(0.25)
                lastIncrementSize = False
                
        
        stepSuccess = stepProgress >= (1-1e-15)
        finishedTime = stepProgress * stepLength
        
        self.journal.message("Time in elements:       {:} s".format(computationTimeInElements), self.identification, level=1)
        self.journal.message("Time in matrix const.  :{:} s".format(computationTimeMatrixConstruction), self.identification, level=1)
        self.journal.message("Time in linear solver:  {:} s".format(computationTimeInEqSystem), self.identification, level=1)
            
        return stepSuccess, U, P, finishedTime
    
    def computeElements(self, U, dU, time, dT, pNewDT, P, V, I, J):
        """ Loop over all elements, and evalute them. 
        Note that ABAQUS style is employed: element(Un+1, dUn+1) 
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration """
            
        P[:] = 0.0
        UN1 = dU + U
        
        for el in self.elements.values():
            idxInVIJ = self.elementToIndexInVIJMap[el]
            Ke = V[idxInVIJ : idxInVIJ+el.sizeKe]
            Pe = np.zeros(el.nDofPerEl)
            idcsInPUdU = I[idxInVIJ : idxInVIJ+el.nDofPerEl]
            
            el.computeYourself(Ke, 
                               Pe, 
                               UN1[ idcsInPUdU ], 
                               dU [ idcsInPUdU ], 
                               time, dT, pNewDT)
            if pNewDT[0] <= 1.0:
                break 
            
            # global effort vector is assembled directly
            P[ idcsInPUdU ] += Pe
            
        return P, V, pNewDT
    
    def applyDirichletK(self, K, dirichlet):
        """ Apply the dirichlet bcs on the global stiffnes matrix
        -> is called by solveStep() before solving the global sys.
        http://stackoverflow.com/questions/12129948/scipy-sparse-set-row-to-zeros"""
            
        for row in dirichlet['indices']:
            K.data[K.indptr[row]:K.indptr[row+1]] = 0.0
            K[row, row] = 1.0
        K.eliminate_zeros()
        return K
              
    def applyDirichlet(self, stepProgress, R, dirichlet):
        """ Apply the dirichlet bcs on the Residual vector
        -> is called by solveStep() before solving the global sys."""
            
        indices = dirichlet['indices']
        R[indices] = dirichlet['delta'] * stepProgress
        return R
    
    def checkConvergency(self, R, ddU, iterationCounter):
        """ Check the convergency, indivudually for each field,
        similar to ABAQUS based on the current total residual and the flow correction
        -> is called by solveStep() to decice wether to continue iterating or stop"""
        
        iterationMessage = ''
        convergedAtAll  = True
        
        if iterationCounter < 15: # standard tolerance set
            effortResidualTolerances = self.effortResidualTolerances
        else: # alternative tolerance set
            effortResidualTolerances = self.effortResidualTolerancesAlt
        
        for field, fieldIndices in self.fieldIndices.items():
            effortResidual =    np.linalg.norm(R[fieldIndices] , np.inf)
            flowCorrection =    np.linalg.norm(ddU[fieldIndices] , np.inf) if ddU is not None else 0.0
            convergedEffort =   True if (effortResidual < effortResidualTolerances[field])  else False
            convergedFlow =     True if (flowCorrection < self.flowCorrectionTolerances[field])  else False
                                     
            iterationMessage += self.iterationMessageTemplate.format(
                                 effortResidual, 
                                 '✓' if convergedEffort  else ' ',
                                 flowCorrection,
                                 '✓' if convergedFlow else ' ',
                                 )
            # converged if residual and flowCorrection are smaller than tolerance
            convergedAtAll = convergedAtAll and convergedFlow and convergedEffort
            
        self.journal.message(iterationMessage, self.identification)     
        return convergedAtAll
    
    def generateVIJ(self, elements):
        """ Initializes the V vector and generates I, J entries for each element,
        based on i) its (global) nodes ii) its dofLayout. Furthermore, 
        a dictionary with the mapping of each element to its index in VIJ 
        is created.
        -> is called by __init__() """
        
        V = np.zeros(self.sizeVIJ)
        I = np.zeros_like(V, dtype=np.int)
        J = np.zeros_like(V, dtype=np.int)
        idxInVIJ = 0
        elementToIndexInVIJMap = {}
        for el in elements.values():
            destList = np.asarray([i for iNode, node in enumerate(el.nodes) 
                                        for nodeField in el.fields[iNode] 
                                            for i in node.fields[nodeField]]) 
    
            elementToIndexInVIJMap[el] = idxInVIJ        
                                  
            # looks like black magic, but it's an efficient way to generate all indices of Ke in K:
            elDofLocations = np.tile(destList[ el.dofIndicesPermutation  ], (destList.shape[0], 1) )
            I[idxInVIJ : idxInVIJ+el.sizeKe] = elDofLocations.ravel()
            J[idxInVIJ : idxInVIJ+el.sizeKe] = elDofLocations.ravel('F')
            idxInVIJ += el.sizeKe
            
        return V, I, J, elementToIndexInVIJMap
             
