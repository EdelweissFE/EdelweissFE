#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 20:37:35 2017

@author: matthias

Standard nonlinear, implicit static solver.

"""
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from fe.utils.incrementgenerator import IncrementGenerator
from fe.utils.exceptions import ReachedMaxIncrements, ReachedMaxIterations, ReachedMinIncrementSize, CutbackRequest, DivergingSolution
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
    defaultMaxGrowingIter = 3
    
    def __init__(self, jobInfo, modelInfo, journal, fieldOutputController, outputmanagers):
        self.nodes =        modelInfo['nodes']
        self.elements =     modelInfo['elements']
        self.nodeSets =     modelInfo['nodeSets']
        self.elementSets =  modelInfo['elementSets']
        self.constraints =  modelInfo['constraints']
        self.fieldIndices = jobInfo['fieldIndices']
        self.residualHistories = dict.fromkeys( self.fieldIndices )
        
        self.fieldCorrectionTolerances =    jobInfo['fieldCorrectionTolerance']
        self.fluxResidualTolerances =       jobInfo['fluxResidualTolerance']
        self.fluxResidualTolerancesAlt =    jobInfo['fluxResidualToleranceAlternative']
        
        # create headers for formatted output of solver
        nFields = len(self.fieldIndices.keys())
        self.iterationHeader = ("{:^25}"*nFields).format(*self.fieldIndices.keys())
        self.iterationHeader2 = (" {:<10}  {:<10}  ").format('||R||∞','||ddU||∞') *nFields
        self.iterationMessageTemplate = "{:11.2e}{:1}{:11.2e}{:1} "
        
        self.nDof = jobInfo['numberOfDofs']
        self.journal = journal
        self.fieldOutputController = fieldOutputController
        self.outputmanagers = outputmanagers 
        
        self.sizeVIJ = 0
        self.sizeNDofElementWise = 0
        
        for el in self.elements.values():
            self.sizeVIJ += el.sizeKe
        
        for constraint in self.constraints.values():
            self.sizeVIJ += constraint.sizeStiffness

        # create indices map to elements; V, I, J are of type np vectors
        # elementToIndexInVIJMap is a dictionary of {element : index in VIJ vectors} 
        V, I, J, elementToIndexInVIJMap , constraintToIndexInVIJMap = self.generateVIJ(self.elements, self.constraints)
        self.V = V
        self.I = I
        self.J = J
        self.elementToIndexInVIJMap = elementToIndexInVIJMap # element  -> V[ .... idx ..  ]
        self.constraintToIndexInVIJMap = constraintToIndexInVIJMap
        
        # for the Abaqus like convergence test, the number of dofs 'element-wise' is needed:
        # = Σ_elements Σ_nodes ( nDof (field) )
        self.fieldNDofElementWise = {}
        for field in self.fieldIndices.keys():
             self.fieldNDofElementWise[field] = np.sum(
                     [ len(node.fields[field]) for el in self.elements.values() 
                                                 for node in (el.nodes) if field in node.fields])
        
    def initialize(self):
        """ Initialize the solver and return the 2 vectors for field (U) and flux (P) """
        
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
        
        R =             np.copy(P)          # global Residual vector for the Newton iteration
        F =             np.zeros_like(P)    # accumulated Flux vector 
        Pdeadloads =    np.zeros_like(P)
        
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
        maxGrowingIter = step.get('maxGrowIter', self.defaultMaxGrowingIter)
        
        dU = np.zeros(numberOfDofs)
        
        dirichlets =        list(stepActions['dirichlet'].values()) + list(stepActions['shotcreteshellmaster'].values())
        distributedLoads =  stepActions['distributedload'].values()
        concentratedLoads = stepActions['nodeforces'].values()
        
        geostatics = stepActions['geostatic'].values()
        activeGeostatics = [g for g in geostatics if g.active]
        isGeostaticStep = any([g.active for g in geostatics] )
        linearConstraints = [c for c in self.constraints.values() if c.linearConstraint ]
        
        for constraint in linearConstraints:
            # linear constraints are independent on the solution
            idxInVIJ = self.constraintToIndexInVIJMap[constraint]
            V[idxInVIJ : idxInVIJ+constraint.sizeStiffness] = constraint.K.ravel()
            
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
                incrementResidualHistory = dict.fromkeys( self.fieldIndices, (0.0, 0 ) )
                spatialAveragedFluxes = dict.fromkeys(self.fieldIndices, 0.0)
                
                stepTimes = np.array([stepTime, totalTime])
                
                if extrapolation == 'linear' and lastIncrementSize:
                    dU *= (incrementSize/lastIncrementSize) 
                    dU = self.applyDirichlet(increment, dU, dirichlets)      
                    extrapolatedIncrement = True
                else:
                    extrapolatedIncrement = False
                    dU[:] = 0.0    
                    
                # Deadloads: only once per increment
                if concentratedLoads or distributedLoads:
                    Pdeadloads[:] = 0.0
                    for cLoad in concentratedLoads: Pdeadloads = cLoad.applyOnP(Pdeadloads, increment) 
                    Pdeadloads = self.computeDistributedLoads(distributedLoads, Pdeadloads, 
                                                              I, stepTimes, dT, increment)
    
                try:
                    while True:
                        for geostatic in activeGeostatics: geostatic.apply() 
                        
                        P, V, F = self.computeElements(U, dU, stepTimes, dT, P, V, I, J, F)
                        
                        if concentratedLoads or distributedLoads:
                            P += Pdeadloads
                            
                        R[:] = P
                        
                        if iterationCounter == 0 and not extrapolatedIncrement and dirichlets :
                            # first iteration? apply dirichlet bcs and unconditionally solve
                            R = self.applyDirichlet(increment, R, dirichlets)
                        else:
                            # iteration cycle 1 or higher, time to check the convergency
                            for dirichlet in dirichlets: 
                                R[dirichlet.indices] = 0.0 
                            for constraint in self.constraints.values(): 
                                R[constraint.globalDofIndices] = 0.0 # currently no external loads on rbs possible
                            
                            for field, nDof in self.fieldNDofElementWise.items():
                                spatialAveragedFluxes[field] = np.linalg.norm(F[self.fieldIndices[field]],1) / nDof 
                            
                            if self.checkConvergency(R, ddU, spatialAveragedFluxes, iterationCounter, incrementResidualHistory) :
                                break
                            
                            for previousFluxResidual, nGrew in incrementResidualHistory.values():
                                if nGrew > maxGrowingIter: 
                                    raise DivergingSolution('Residual grew {:} times, cutting back'.format(maxGrowingIter))
                            
                        if iterationCounter == maxIter:
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
                    
                except (ReachedMaxIterations, 
                        DivergingSolution) as e:
                    self.journal.message(str(e), self.identification, 1)
                    incGen.discardAndChangeIncrement(0.25)
                    lastIncrementSize = False
                    
                else: 
                    U += dU
                    lastIncrementSize = incrementSize
                    if iterationCounter >= criticalIter: 
                        incGen.preventIncrementIncrease()
                        
                    for el in self.elements.values(): 
                        el.acceptLastState()
                        
                    self.journal.message("Converged in {:} iteration(s)".format(iterationCounter), self.identification, 1) 
                    
                    self.fieldOutputController.finalizeIncrement(U, P, increment)
                    for man in self.outputmanagers:
                        man.finalizeIncrement(U, P, increment)
                        
        except (ReachedMaxIncrements, ReachedMinIncrementSize) as e:
            self.journal.errorMessage("Incrementation failed", self.identification)
            
        except KeyboardInterrupt:
            print('')
            self.journal.message("Interrupted by user", self.identification)
        
        else:
            if isGeostaticStep: U = self.resetDisplacements(U)  # reset all displacements, if the present step is a geostatic step
            
            for stepActionType in stepActions.values():
                for action in stepActionType.values():
                    action.finishStep()
                    
        finally:
            finishedTime = time + stepProgress * stepLength
            self.journal.printTable([ ("Time in {:}".format(k), " {:10.4f}s".format(v)) for k, v in self.computationTimes.items() ], self.identification)
            
        return ((1.0 - stepProgress) < 1e-14) , U, P, finishedTime
    
    def computeElements(self, U, dU, time, dT, P, V, I, J, F):
        """ Loop over all elements, and evalute them. 
        Note that ABAQUS style is employed: element(Un+1, dUn+1) 
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration """
        
        tic = getCurrentTime()
        P[:] =      0.0
        F[:] = 0.0
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
            
            # global force vector is assembled directly
            P[ idcsInPUdU ] += Pe
            F[ idcsInPUdU ] += np.abs(Pe)
        
        toc = getCurrentTime()
        self.computationTimes['elements'] += toc - tic
        return P, V, F
    
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
        http://stackoverflux.com/questions/12129948/scipy-sparse-set-row-to-zeros"""
            
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
            R[dirichlet.indices] = dirichlet.getDelta(increment)
        
        toc = getCurrentTime()
        self.computationTimes['dirichlet R'] += toc - tic
        
        return R
    
    def checkConvergency(self, R, ddU, spatialAveragedFluxes, iterationCounter, residualHistory):
        """ Check the convergency, indivudually for each field,
        similar to ABAQUS based on the current total residual and the flux correction
        -> is called by solveStep() to decide wether to continue iterating or stop"""
        
        tic = getCurrentTime()
        
        iterationMessage = ''
        convergedAtAll  = True
        
        if iterationCounter < 15: # standard tolerance set
            fluxResidualTolerances = self.fluxResidualTolerances
        else: # alternative tolerance set
            fluxResidualTolerances = self.fluxResidualTolerancesAlt
        
        for field, fieldIndices in self.fieldIndices.items():
            fluxResidual =       np.linalg.norm(R[fieldIndices] , np.inf)
            fieldCorrection =    np.linalg.norm(ddU[fieldIndices] , np.inf) if ddU is not None else 0.0
            convergedCorrection = True if (fieldCorrection < self.fieldCorrectionTolerances[field]) else False
            convergedFlux =       True if (fluxResidual <= fluxResidualTolerances[field] * spatialAveragedFluxes[field]  )  else False
            
            previousFluxResidual, nGrew = residualHistory[field]
            if fluxResidual > previousFluxResidual:
                nGrew += 1
            residualHistory[field]  = (fluxResidual, nGrew)
                                     
            iterationMessage += self.iterationMessageTemplate.format(
                                 fluxResidual, 
                                 '✓' if convergedFlux  else ' ',
                                 fieldCorrection,
                                 '✓' if convergedCorrection else ' ',
                                 )
            # converged if residual and fieldCorrection are smaller than tolerance
            convergedAtAll = convergedAtAll and convergedCorrection and convergedFlux
            
            
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
    
    def tocsr(self, coo, copy=False):
        from scipy.sparse.sputils import  get_index_dtype, upcast
        from scipy.sparse._sparsetools import coo_tocsr
        M,N = coo.shape
        idx_dtype = get_index_dtype((coo.row, coo.col),
                                    maxval=max(coo.nnz, N))
        row = coo.row.astype(idx_dtype, copy=False)
        col = coo.col.astype(idx_dtype, copy=False)

        indptr = np.empty(M + 1, dtype=idx_dtype)
        indices = np.empty_like(col, dtype=idx_dtype)
        data = np.empty_like(coo.data, dtype=upcast(coo.dtype))

        coo_tocsr(M, N, coo.nnz, row, col, coo.data,
                  indptr, indices, data)

        x = csr_matrix((data, indices, indptr), shape=coo.shape)
        if not coo.has_canonical_format:
            x.sum_duplicates()
        return x
    
    def assembleStiffness(self, V, I, J, shape):
        """ Construct a CSR matrix from VIJ """
        tic =  getCurrentTime()
#        K  = coo_matrix( (V, (I,J)), shape).tocsr()
        K = self.tocsr(coo_matrix( (V, (I,J)), shape))
        toc =  getCurrentTime()
        self.computationTimes['CSR generation'] += toc - tic
        return K
    
    def generateVIJ(self, elements, constraints):
        """ Initializes the V vector and generates I, J entries for each element,
        based on i) its (global) nodes ii) its dofLayout. Furthermore, 
        a dictionary with the mapping of each element to its index in VIJ 
        is created.
        The same is done for each constraint
        -> is called by __init__() """
        
        V = np.zeros(self.sizeVIJ)
        I = np.zeros_like(V, dtype=np.int)
        J = np.zeros_like(V, dtype=np.int)
        idxInVIJ = 0
        elementToIndexInVIJMap = {}
        
        for el in elements.values():
            destList = np.asarray([i for iNode, node in enumerate(el.nodes) # for each node of the element..
                                        for nodeField in el.fields[iNode]  # for each field of this node
                                            for i in node.fields[nodeField]])  # the index in the global system
    
    
            elementToIndexInVIJMap[el] = idxInVIJ        
                                  
            # looks like black magic, but it's an efficient way to generate all indices of Ke in K:
            elDofLocations = np.tile(destList[ el.dofIndicesPermutation  ], (destList.shape[0], 1) )
            I[idxInVIJ : idxInVIJ+el.sizeKe] = elDofLocations.ravel()
            J[idxInVIJ : idxInVIJ+el.sizeKe] = elDofLocations.ravel('F')
            idxInVIJ += el.sizeKe
            
        constraintToIndexInVIJMap = {}
        for constraint in constraints.values():
            destList = constraint.globalDofIndices
            constraintToIndexInVIJMap[constraint] = idxInVIJ        
            constraintDofLocations = np.tile( destList, (destList.shape[0], 1) )
            I[idxInVIJ : idxInVIJ + constraint.sizeStiffness] = constraintDofLocations.ravel()
            J[idxInVIJ : idxInVIJ + constraint.sizeStiffness] = constraintDofLocations.ravel('F')
            idxInVIJ += constraint.sizeStiffness
        
        return V, I, J, elementToIndexInVIJMap, constraintToIndexInVIJMap
    
    def resetDisplacements(self, U):
        """ -> May be called at the end of a geostatic step, the zero all displacments after
        applying the geostatic stress """
        
        self.journal.printSeperationLine()
        self.journal.message("Geostatic step, resetting displacements", self.identification, level=1)
        self.journal.printSeperationLine()
        U[ self.fieldIndices['displacement'] ] = 0.0 
        return U
    
