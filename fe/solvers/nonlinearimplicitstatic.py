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
from fe.utils.exceptions import ReachedMaxIncrements, ReachedMaxIterations, ReachedMinIncrementSize, CutbackRequest
from time import time as getCurrentTime
from collections import defaultdict

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
        
    def solveStep(self, step, time, stepActions, stepOptions, U, P):
        """ Public interface to solve for an ABAQUS like step
        returns: boolean Success, U vector, P vector, and the new current total time """
    
        self.computationTimes = defaultdict(lambda : 0.0)
            
        V = self.V
        J = self.J
        I = self.I
        
        R = np.copy(P)
        
        numberOfDofs = self.nDof
        stepLength = step.get('stepLength', 1.0)
        extrapolation = stepOptions['NISTSolver'].get('extrapolation', 'linear')

        incGen = IncrementGenerator(time, 
                                    stepLength, 
                                    step.get('maxInc', self.defaultMaxInc), 
                                    step.get('minInc', self.defaultMinInc), 
                                    step.get('maxNumInc', self.defaultMaxNumInc), 
                                    self.journal)
        
        maxIter = step.get('maxIter', self.defaultMaxIter)
        criticalIter = step.get('crititcalIter', self.defaultCriticalIter)
        
        dU = np.zeros(numberOfDofs)
        
        dirichlets = stepActions['dirichlet'].values()
        distributedLoads = stepActions['distributedload'].values()
        concentratedLoads = stepActions['nodeforces'].values()
        geostatics = stepActions['geostatic'].values()
        isGeostaticStep = any([g.active for g in geostatics] )
        
        # initialize pNewDt as np.array, for bypassing to external c++ code (UEL)
        lastIncrementSize = False
        
        try:
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
                
                ddU = None
                iterationCounter = 0
                stepTimes = np.array([stepTime, totalTime])
                
                if extrapolation == 'linear' and lastIncrementSize:
                    dU *= (incrementSize/lastIncrementSize) 
                    dU = self.applyDirichlet(increment, dU, dirichlets)      
                    extrapolatedIncrement = True
                else:
                    extrapolatedIncrement = False
                    dU[:] = 0.0        
    
                try:
                    while True:
                        
                        for geostatic in geostatics:
                            geostatic.apply()
                        
                        P, V = self.computeElements(U, dU, stepTimes, dT, P, V, I, J,)
                            
                        for cLoad in concentratedLoads: 
                            P = cLoad.applyOnP(P, increment) 
                        self.computeDistributedLoads(distributedLoads, P, I, stepTimes, dT, increment)   
                        
                        R[:] = P
                        
                        if iterationCounter == 0 and not extrapolatedIncrement and dirichlets :
                            # first iteration? apply dirichlet bcs and unconditionally solve
                            R = self.applyDirichlet(increment, R, dirichlets)
                        else:
                            # iteration cycle 1 or higher, time to check the convergency
                            for dirichlet in dirichlets: R[dirichlet.indices] = 0.0 # only entries not affected by dirichlet bcs contribute to the residual
                            
                            # check convergency 
                            if self.checkConvergency(R, ddU, iterationCounter):
                                break
                            
                        if iterationCounter == maxIter-1:
                            raise  ReachedMaxIterations("Reached max. iterations in current increment, cutting back")
                        
                        K = self.assembleStiffness(V, I, J, shape=(numberOfDofs, numberOfDofs) )
                        K = self.applyDirichletK(K, dirichlets)
                        ddU = self.linearSolve(K, R, )
                        dU += ddU
                        iterationCounter += 1
                    
                except CutbackRequest as e:
                    self.journal.message(str(e), self.identification, 1)
                    incGen.discardAndChangeIncrement( e.cutbackSize if e.cutbackSize > 0.25 else 0.25 )
                    lastIncrementSize = False
                    
                except ReachedMaxIterations as e:
                    self.journal.message(str(e), self.identification, 1)
                    incGen.discardAndChangeIncrement(0.25)
                    lastIncrementSize = False
                    
                else: 
                    #converged
                    U += dU
                    lastIncrementSize = incrementSize
                    if iterationCounter >= criticalIter: 
                        incGen.preventIncrementIncrease()
                        
                    for el in self.elements.values(): 
                        el.acceptLastState()
                        
                    self.journal.message("Converged in {:} iteration(s)".format(iterationCounter), self.identification, 1) 
                    
                    for man in self.outputmanagers:
                        man.finalizeIncrement(U, P, increment)
                        
        except (ReachedMaxIncrements, ReachedMinIncrementSize) as e:
            self.journal.errorMessage("Incrementation failed", self.identification)
            
        except KeyboardInterrupt:
            print('')
            self.journal.message("Interrupted by user", self.identification)
            
        else:
            # reset all displacements, if the present step is a geostatic step
            if isGeostaticStep: U = self.resetDisplacements(U)
                
            for stepActionType in stepActions.values():
                for action in stepActionType.values():
                    action.finishStep()
                    
        finally:
            
            finishedTime = time + stepProgress * stepLength
            for section, time in self.computationTimes.items():
                self.journal.message("Time in {:<30}: {:}s".format(section, time), self.identification, level=1)
            
        return 1.0 - stepProgress < 1e-15 , U, P, finishedTime
    
    def computeElements(self, U, dU, time, dT, P, V, I, J):
        """ Loop over all elements, and evalute them. 
        Note that ABAQUS style is employed: element(Un+1, dUn+1) 
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration """
        tic = getCurrentTime()
        
        P[:] = 0.0
        UN1 = dU + U
        pNewDT = np.array([1e36])
        
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
                raise CutbackRequest("An element requests for a cutback", pNewDT[0])
            
            # global effort vector is assembled directly
            P[ idcsInPUdU ] += Pe
        
        toc = getCurrentTime()
        self.computationTimes['elements'] += toc - tic
        return P, V
    
    def computeDistributedLoads(self, distributedLoads, P, I, time, dT, increment):
        """ Loop over all elements, and evalute them. 
        Note that ABAQUS style is employed: element(Un+1, dUn+1) 
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration """
        
        tic = getCurrentTime()
        for dLoad in distributedLoads:
            magnitude = dLoad.getCurrentMagnitude(increment)
            for faceID, elements in dLoad.surface.items():
                for el in elements:
                    idxInVIJ = self.elementToIndexInVIJMap[el]
                    Pe = np.zeros(el.nDofPerEl)
                    idcsInPUdU = I[idxInVIJ : idxInVIJ+el.nDofPerEl]
                    Pe[:] = 0.0 
                    el.computeDistributedLoad(dLoad.loadType,
                                              Pe,
                                              faceID,
                                              magnitude, 
                                              time, dT)
                    P[idcsInPUdU] += Pe                    
            
        toc = getCurrentTime()
        self.computationTimes['distributed loads'] += toc - tic
        return P
    
    def applyDirichletK(self, K, dirichlets):
        """ Apply the dirichlet bcs on the global stiffnes matrix
        -> is called by solveStep() before solving the global sys.
        http://stackoverflow.com/questions/12129948/scipy-sparse-set-row-to-zeros"""
            
        tic = getCurrentTime()
        for dirichlet in dirichlets:
            for row in dirichlet.indices:
                K.data[K.indptr[row]:K.indptr[row+1]] = 0.0
                K[row, row] = 1.0
        K.eliminate_zeros()
        
        toc = getCurrentTime()
        self.computationTimes['dirichlet K'] += toc - tic
        
        return K
              
    def applyDirichlet(self, increment, R, dirichlets):
        """ Apply the dirichlet bcs on the Residual vector
        -> is called by solveStep() before solving the global sys."""
        tic = getCurrentTime()
        for dirichlet in dirichlets:
            indices = dirichlet.indices#['indices']
            R[indices] = dirichlet.getDelta(increment)

        toc = getCurrentTime()
        self.computationTimes['dirichlet R'] += toc - tic
        
        return R
    
    def checkConvergency(self, R, ddU, iterationCounter):
        """ Check the convergency, indivudually for each field,
        similar to ABAQUS based on the current total residual and the flow correction
        -> is called by solveStep() to decide wether to continue iterating or stop"""
        
        tic = getCurrentTime()
        
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
        
        toc = getCurrentTime()
        self.computationTimes['convergence check'] += toc - tic
        
        return convergedAtAll
    
    def linearSolve(self, A, b):
        """ Interface to the linear eq. system solver """
        
        tic =  getCurrentTime()
        ddU = spsolve(A, b )
        toc =  getCurrentTime()
        self.computationTimes['linear solve'] += toc - tic
        return ddU
    
    def assembleStiffness(self, V, I, J, shape):
        """ Construct a CSR matrix from VIJ """
        
        tic =  getCurrentTime()
        K = coo_matrix( (V, (I,J)), shape).tocsr()
        toc =  getCurrentTime()
        self.computationTimes['CSR generation'] += toc - tic
        return K
    
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
    
    def resetDisplacements(self, U):
        """ -> May be called at the end of a geostatic step, the zero all displacments after
        applying the geostatic stress """
        
        self.journal.printSeperationLine()
        self.journal.message("Geostatic step, resetting displacements", self.identification, level=1)
        self.journal.printSeperationLine()
        U[ self.fieldIndices['displacement'] ] = 0.0 
        return U
    