#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  ---------------------------------------------------------------------
#
#  _____    _      _              _         _____ _____
# | ____|__| | ___| |_      _____(_)___ ___|  ___| ____|
# |  _| / _` |/ _ \ \ \ /\ / / _ \ / __/ __| |_  |  _|
# | |__| (_| |  __/ |\ V  V /  __/ \__ \__ \  _| | |___
# |_____\__,_|\___|_| \_/\_/ \___|_|___/___/_|   |_____|
#
#
#  Unit of Strength of Materials and Structural Analysis
#  University of Innsbruck,
#  2017 - today
#
#  Matthias Neuner matthias.neuner@uibk.ac.at
#
#  This file is part of EdelweissFE.
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  The full text of the license can be found in the file LICENSE.md at
#  the top level directory of EdelweissFE.
#  ---------------------------------------------------------------------
# Created on Sun Jan  8 20:37:35 2017

# @author: Matthias Neuner

import numpy as np
from fe.utils.incrementgenerator import IncrementGenerator
from fe.utils.exceptions import (
    ReachedMaxIncrements,
    ReachedMaxIterations,
    ReachedMinIncrementSize,
    CutbackRequest,
    DivergingSolution,
    ConditionalStop,
)
from time import time as getCurrentTime
from collections import defaultdict
from fe.config.linsolve import getLinSolverByName, getDefaultLinSolver
from fe.utils.csrgenerator import CSRGenerator
from fe.stepactions.base.stepactionbase import StepActionBase
from fe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from fe.utils.dofmanager import DofVector, VIJSystemMatrix
from fe.utils.fieldoutput import FieldOutputController
from fe.constraints.base.constraintbase import ConstraintBase
from scipy.sparse import csr_matrix
from numpy import ndarray


class NIST:
    """This is the Nonlinear Implicit STatic -- solver.

    Parameters
    ----------
    jobInfo
        A dictionary containing the job information.
    journal
        The journal instance for logging.
    """

    identification = "NISTSolver"

    defaultMaxInc = 1.0
    defaultMinInc = 1e-4
    defaultMaxNumInc = 1000
    defaultMaxIter = 10
    defaultCriticalIter = 5
    defaultMaxGrowingIter = 10

    def __init__(self, jobInfo, journal):

        self.linSolver = getLinSolverByName(jobInfo["linsolver"]) if "linsolver" in jobInfo else getDefaultLinSolver()

        self.theDofManager = jobInfo["dofManager"]
        self.fieldCorrectionTolerances = jobInfo["fieldCorrectionTolerance"]
        self.fluxResidualTolerances = jobInfo["fluxResidualTolerance"]
        self.fluxResidualTolerancesAlt = jobInfo["fluxResidualToleranceAlternative"]

        # create headers for formatted output of solver
        nFields = len(self.theDofManager.indicesOfFieldsInDofVector.keys())
        self.iterationHeader = ("{:^25}" * nFields).format(*self.theDofManager.indicesOfFieldsInDofVector.keys())
        self.iterationHeader2 = (" {:<10}  {:<10}  ").format("||R||∞", "||ddU||∞") * nFields
        self.iterationMessageTemplate = "{:11.2e}{:1}{:11.2e}{:1} "

        self.journal = journal

        self.systemMatrix = self.theDofManager.constructVIJSystemMatrix()

        self.csrGenerator = CSRGenerator(self.systemMatrix)

        self.residualHistories = dict.fromkeys(self.theDofManager.indicesOfFieldsInDofVector)

        self.extrapolation = "linear"

    def initialize(self):
        """Initialize the solver and return the 2 vectors for field (U) and flux (P)."""

        U = self.theDofManager.constructDofVector()
        P = self.theDofManager.constructDofVector()

        return U, P

    def solveStep(
        self,
        stepNumber: int,
        step: dict,
        time: float,
        stepActions: dict[str, StepActionBase],
        modelInfo: dict,
        U: DofVector,
        P: DofVector,
        fieldOutputController: FieldOutputController,
        outputmanagers: dict[str, OutputManagerBase],
    ) -> tuple[bool, DofVector, DofVector, float]:
        """Public interface to solve for a step.

        Parameters
        ----------
        stepNumber
            The step number.
        step
            The dictionary containing the step definition.
        time
            The time at beginning of the step.
        stepActions
            The dictionary containing all step actions.
        modelInfo
            A dictionary containing the model tree.
        U
            The current solution vector.
        P
            The current reaction vector.
        fieldOutputController
            The field output controller.

        Returns
        -------
        tuple[bool,DofVector,DofVector,float]
            - Truth value of success.
            - The new solution vector.
            - The new reaction vector.
            - The new current time.

        """

        self.computationTimes = defaultdict(lambda: 0.0)

        K = self.systemMatrix

        stepLength = step.get("stepLength", 1.0)

        extrapolation = self.extrapolation
        try:
            extrapolation = stepActions["options"]["NISTSolver"]["extapolation"]
        except KeyError:
            pass

        incGen = IncrementGenerator(
            time,
            stepLength,
            step.get("maxInc", self.defaultMaxInc),
            step.get("minInc", self.defaultMinInc),
            step.get("maxNumInc", self.defaultMaxNumInc),
            self.journal,
        )

        maxIter = step.get("maxIter", self.defaultMaxIter)
        criticalIter = step.get("criticalIter", self.defaultCriticalIter)
        maxGrowingIter = step.get("maxGrowIter", self.defaultMaxGrowingIter)

        dU = self.theDofManager.constructDofVector()

        elements = modelInfo["elements"]

        constraints = modelInfo["constraints"]

        lastIncrementSize = False
        success = False

        self.applyStepActionsAtStepStart(U, P, stepActions)

        try:
            for increment in incGen.generateIncrement():
                incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment

                statusInfoDict = {
                    "step": stepNumber,
                    "inc": incNumber,
                    "iters": None,
                    "converged": False,
                    "time inc": dT,
                    "time": totalTime + dT,
                    "notes": "",
                }

                self.journal.printSeperationLine()
                self.journal.message(
                    "increment {:}: {:8f}, {:8f}; time {:10f} to {:10f}".format(
                        incNumber, incrementSize, stepProgress, totalTime, totalTime + dT
                    ),
                    self.identification,
                    level=1,
                )
                self.journal.message(self.iterationHeader, self.identification, level=2)
                self.journal.message(self.iterationHeader2, self.identification, level=2)

                try:
                    U, dU, P, iterationCounter, incrementResidualHistory = self.solveIncrement(
                        U,
                        dU,
                        P,
                        K,
                        elements,
                        stepActions,
                        constraints,
                        increment,
                        lastIncrementSize,
                        extrapolation,
                        maxIter,
                        maxGrowingIter,
                    )

                except CutbackRequest as e:
                    self.journal.message(str(e), self.identification, 1)
                    incGen.discardAndChangeIncrement(max(e.cutbackSize, 0.25))
                    lastIncrementSize = False

                    statusInfoDict["iters"] = np.inf
                    statusInfoDict["notes"] = str(e)

                    for man in outputmanagers:
                        man.finalizeFailedIncrement(statusInfoDict=statusInfoDict)

                except (ReachedMaxIterations, DivergingSolution) as e:
                    self.journal.message(str(e), self.identification, 1)
                    incGen.discardAndChangeIncrement(0.25)
                    lastIncrementSize = False

                    statusInfoDict["iters"] = np.inf
                    statusInfoDict["notes"] = str(e)

                    for man in outputmanagers:
                        man.finalizeFailedIncrement(statusInfoDict=statusInfoDict)

                else:
                    lastIncrementSize = incrementSize
                    if iterationCounter >= criticalIter:
                        incGen.preventIncrementIncrease()

                    for el in elements.values():
                        el.acceptLastState()

                    self.journal.message(
                        "Converged in {:} iteration(s)".format(iterationCounter), self.identification, 1
                    )

                    statusInfoDict["iters"] = iterationCounter
                    statusInfoDict["converged"] = True

                    fieldOutputController.finalizeIncrement(U, P, increment)
                    for man in outputmanagers:
                        man.finalizeIncrement(U, P, increment, statusInfoDict=statusInfoDict)

        except (ReachedMaxIncrements, ReachedMinIncrementSize):
            success = False
            self.journal.errorMessage("Incrementation failed", self.identification)

        except KeyboardInterrupt:
            success = False
            print("")
            self.journal.message("Interrupted by user", self.identification)

        except ConditionalStop:
            success = True
            self.journal.message("Conditional Stop", self.identification)
            self.applyStepActionsAtStepEnd(U, P, stepActions)

        else:
            success = True
            self.applyStepActionsAtStepEnd(U, P, stepActions)

        finally:
            finishedTime = time + stepProgress * stepLength
            self.journal.printTable(
                [("Time in {:}".format(k), " {:10.4f}s".format(v)) for k, v in self.computationTimes.items()],
                self.identification,
            )

        return success, U, P, finishedTime

    def solveIncrement(
        self,
        U_n: DofVector,
        dU: DofVector,
        P: DofVector,
        K: VIJSystemMatrix,
        elements: dict,
        stepActions: list,
        constraints: dict,
        increment: tuple,
        lastIncrementSize: float,
        extrapolation: str,
        maxIter: int,
        maxGrowingIter: int,
    ) -> tuple[DofVector, DofVector, DofVector, int, dict]:
        """Standard Newton-Raphson scheme to solve for an increment.

        Parameters
        ----------
        Un
            The old solution vector.
        dU
            The old solution increment.
        P
            The old reaction vector.
        K
            The system matrix to be used.
        elements
            The dictionary containing all elements.
        stepActions
            The list of active step actions.
        constraints
            The dictionary containing all elements.
        increment
            The increment.
        lastIncrementSize
            The size of the previous increment.
        extrapolation
            The type of extrapolation to be used.
        maxIter
            The maximum number of iterations to be used.
        maxGrowingIter
            The maximum number of growing residuals until the Newton-Raphson is terminated.

        Returns
        -------
        tuple[DofVector,DofVector,DofVector,int,dict]
            A tuple containing
                - the new solution vector
                - the solution increment
                - the new reaction vector
                - the number of required iterations
                - the history of residuals per field
        """

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        iterationCounter = 0
        incrementResidualHistory = dict.fromkeys(self.theDofManager.indicesOfFieldsInDofVector, (0.0, 0))

        R = self.theDofManager.constructDofVector()
        F = self.theDofManager.constructDofVector()
        PExt = self.theDofManager.constructDofVector()
        U_np = self.theDofManager.constructDofVector()
        ddU = None

        dirichlets = stepActions["dirichlet"].values()
        nodeforces = stepActions["nodeforces"].values()
        distributedLoads = stepActions["distributedload"].values()
        bodyForces = stepActions["bodyforce"].values()

        self.applyStepActionsAtIncrementStart(U_n, P, increment, stepActions)

        dU, isExtrapolatedIncrement = self.extrapolateLastIncrement(
            extrapolation, increment, dU, dirichlets, lastIncrementSize
        )

        while True:
            for geostatic in stepActions["geostatic"].values():
                geostatic.applyAtIterationStart()

            U_np[:] = U_n
            U_np += dU

            P[:] = K[:] = F[:] = PExt[:] = 0.0

            P, K, F = self.computeElements(elements, U_np, dU, P, K, F, increment)
            PExt, K = self.assembleLoads(nodeforces, distributedLoads, bodyForces, U_np, PExt, K, increment)
            PExt, K = self.assembleConstraints(constraints, U_np, dU, PExt, K, increment)

            R[:] = P
            R += PExt

            if iterationCounter == 0 and not isExtrapolatedIncrement and dirichlets:
                # first iteration? apply dirichlet bcs and unconditionally solve
                R = self.applyDirichlet(increment, R, dirichlets)
            else:
                # iteration cycle 1 or higher, time to check the convergence
                for dirichlet in dirichlets:
                    R[dirichlet.indices] = 0.0

                converged, nodesWithLargestResidual = self.checkConvergence(
                    R, ddU, F, iterationCounter, incrementResidualHistory
                )

                if converged:
                    break

                if self.checkDivergingSolution(incrementResidualHistory, maxGrowingIter):
                    self.printResidualOutlierNodes(nodesWithLargestResidual)
                    raise DivergingSolution("Residual grew {:} times, cutting back".format(maxGrowingIter))

                if iterationCounter == maxIter:
                    self.printResidualOutlierNodes(nodesWithLargestResidual)
                    raise ReachedMaxIterations("Reached max. iterations in current increment, cutting back")

            K_ = self.assembleStiffnessCSR(K)
            K_ = self.applyDirichletK(K_, dirichlets)

            ddU = self.linearSolve(K_, R)
            dU += ddU
            iterationCounter += 1

        return U_np, dU, P, iterationCounter, incrementResidualHistory

    def computeDistributedLoads(
        self,
        distributedLoads: list[StepActionBase],
        U_np: DofVector,
        PExt: DofVector,
        K: VIJSystemMatrix,
        increment: tuple,
    ) -> tuple[DofVector, VIJSystemMatrix]:
        """Loop over all distributed loads acting on elements, and evaluate them.
        Assembles into the global external load vector and the system matrix.

        Parameters
        ----------
        distributedLoads
            The list of distributed loads.
        U_np
            The current solution vector.
        PExt
            The external load vector to be augmented.
        K
            The system matrix to be augmented.
        increment
            The increment.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix]
            The augmented load vector and system matrix.
        """

        tic = getCurrentTime()
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        time = np.array([stepTime, totalTime])
        for dLoad in distributedLoads:
            magnitude = dLoad.getCurrentMagnitude(increment)
            for faceID, elements in dLoad.surface.items():
                for el in elements:

                    Ke = K[el]
                    Pe = np.zeros(el.nDof)

                    el.computeDistributedLoad(dLoad.loadType, Pe, Ke, faceID, magnitude, U_np[el], time, dT)

                    PExt[el] += Pe

        toc = getCurrentTime()
        self.computationTimes["distributed loads"] += toc - tic
        return PExt, K

    def computeBodyForces(
        self, bodyForces: list[StepActionBase], U_np: DofVector, PExt: DofVector, K: VIJSystemMatrix, increment: tuple
    ) -> tuple[DofVector, VIJSystemMatrix]:
        """Loop over all body forces loads acting on elements, and evaluate them.
        Assembles into the global external load vector and the system matrix.

        Parameters
        ----------
        distributedLoads
            The list of distributed loads.
        U_np
            The current solution vector.
        PExt
            The external load vector to be augmented.
        K
            The system matrix to be augmented.
        increment
            The increment.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix]
            The augmented load vector and system matrix.
        """

        tic = getCurrentTime()
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        time = np.array([stepTime, totalTime])
        for bForce in bodyForces:
            force = bForce.getCurrentBodyForce(increment)
            for el in bForce.elements:

                Pe = np.zeros(el.nDof)
                Ke = K[el]

                el.computeBodyForce(Pe, Ke, force, U_np[el], time, dT)

                PExt[el] += Pe

        toc = getCurrentTime()
        self.computationTimes["body forces"] += toc - tic
        return PExt, K

    def applyDirichletK(self, K: VIJSystemMatrix, dirichlets: list[StepActionBase]) -> VIJSystemMatrix:
        """Apply the dirichlet bcs on the global stiffness matrix
        Is called by solveStep() before solving the global system.
        http://stackoverflux.com/questions/12129948/scipy-sparse-set-row-to-zeros

        Parameters
        ----------
        K
            The system matrix.
        dirichlets
            The list of dirichlet boundary conditions.

        Returns
        -------
        VIJSystemMatrix
            The modified system matrix.
        """

        tic = getCurrentTime()
        for dirichlet in dirichlets:
            for row in dirichlet.indices:
                K.data[K.indptr[row] : K.indptr[row + 1]] = 0.0

        # K[row, row] = 1.0 @ once, faster than within the loop above:
        diag = K.diagonal()
        diag[np.concatenate([d.indices for d in dirichlets])] = 1.0
        K.setdiag(diag)

        K.eliminate_zeros()

        toc = getCurrentTime()
        self.computationTimes["dirichlet K"] += toc - tic

        return K

    def applyDirichlet(self, increment: tuple, R: DofVector, dirichlets: list[StepActionBase]):
        """Apply the dirichlet bcs on the residual vector
        Is called by solveStep() before solving the global equatuon system.

        Parameters
        ----------
        increment
            The increment.
        R
            The residual vector of the global equation system to be modified.
        dirichlets
            The list of dirichlet boundary conditions.

        Returns
        -------
        DofVector
            The modified residual vector.
        """

        tic = getCurrentTime()
        for dirichlet in dirichlets:
            R[dirichlet.indices] = dirichlet.getDelta(increment)

        toc = getCurrentTime()
        self.computationTimes["dirichlet R"] += toc - tic

        return R

    def checkConvergence(
        self, R: DofVector, ddU: DofVector, F: DofVector, iterationCounter: int, residualHistory: dict
    ) -> tuple[bool, dict]:
        """Check the convergence, individually for each field,
        similar to Abaqus based on the current total flux residual and the field correction
        Is called by solveStep() to decide whether to continue iterating or stop.

        Parameters
        ----------
        R
            The current residual.
        ddU
            The current correction increment.
        F
            The accumulated fluxes.
        iterationCounter
            The current iteration number.
        residualHistory
            The previous residuals.

        Returns
        -------
        tuple[bool,dict]
            - True if converged.
            - The residual histories field wise.

        """

        tic = getCurrentTime()

        iterationMessage = ""
        convergedAtAll = True
        nodesWithLargestResidual = {}

        spatialAveragedFluxes = self.computeSpatialAveragedFluxes(F)

        if iterationCounter < 15:  # standard tolerance set
            fluxResidualTolerances = self.fluxResidualTolerances
        else:  # alternative tolerance set
            fluxResidualTolerances = self.fluxResidualTolerancesAlt

        for field, fieldIndices in self.theDofManager.indicesOfFieldsInDofVector.items():
            fieldResidualAbs = np.abs(R[fieldIndices])

            indexOfMax = np.argmax(fieldResidualAbs)
            fluxResidual = fieldResidualAbs[indexOfMax]

            nodesWithLargestResidual[field] = self.theDofManager.getNodeForIndexInDofVector(indexOfMax)

            fieldCorrection = np.linalg.norm(ddU[fieldIndices], np.inf) if ddU is not None else 0.0

            convergedCorrection = fieldCorrection < self.fieldCorrectionTolerances[field]
            convergedFlux = fluxResidual <= max(fluxResidualTolerances[field] * spatialAveragedFluxes[field], 1e-7)

            previousFluxResidual, nGrew = residualHistory[field]
            if fluxResidual > previousFluxResidual:
                nGrew += 1
            residualHistory[field] = (fluxResidual, nGrew)

            iterationMessage += self.iterationMessageTemplate.format(
                fluxResidual,
                "✓" if convergedFlux else " ",
                fieldCorrection,
                "✓" if convergedCorrection else " ",
            )

            # converged if residual and field correction are smaller than tolerance
            convergedAtAll = convergedAtAll and convergedCorrection and convergedFlux

        self.journal.message(iterationMessage, self.identification)

        toc = getCurrentTime()
        self.computationTimes["convergence check"] += toc - tic

        return convergedAtAll, nodesWithLargestResidual

    def linearSolve(self, A: csr_matrix, b: DofVector) -> ndarray:
        """Solve the linear equation system.

        Parameters
        ----------
        A
            The system matrix in compressed spare row format.
        b
            The right hand side.

        Returns
        -------
        ndarray
            The solution 'x'.
        """

        tic = getCurrentTime()
        ddU = self.linSolver(A, b)
        toc = getCurrentTime()
        self.computationTimes["linear solve"] += toc - tic

        if np.isnan(ddU).any():
            raise DivergingSolution("Obtained NaN in linear solve")

        return ddU

    def assembleStiffnessCSR(self, K: VIJSystemMatrix) -> csr_matrix:
        """Construct a CSR matrix from VIJ format.

        Parameters
        ----------
        K
            The system matrix in VIJ format.
        Returns
        -------
        csr_matrix
            The system matrix in compressed sparse row format.
        """
        tic = getCurrentTime()
        KCsr = self.csrGenerator.updateCSR(K)
        toc = getCurrentTime()
        self.computationTimes["CSR generation"] += toc - tic
        return KCsr

    def computeSpatialAveragedFluxes(self, F: DofVector) -> np.float:
        """Compute the spatial averaged flux for every field
        Is usually called by checkConvergence().

        Parameters
        ----------
        F
            The accumulated flux vector.

        Returns
        -------
        dict[str,np.float]
            A dictioary containg the spatial average fluxes for every field.
        """
        spatialAveragedFluxes = dict.fromkeys(self.theDofManager.indicesOfFieldsInDofVector, 0.0)
        for field, nDof in self.theDofManager.nAccumulatedNodalFluxesFieldwise.items():
            spatialAveragedFluxes[field] = max(
                1e-10, np.linalg.norm(F[self.theDofManager.indicesOfFieldsInDofVector[field]], 1) / nDof
            )

        return spatialAveragedFluxes

    def computeElements(
        self,
        elements: list,
        U_np: DofVector,
        dU: DofVector,
        P: DofVector,
        K: VIJSystemMatrix,
        F: DofVector,
        increment: tuple,
    ) -> tuple[DofVector, VIJSystemMatrix, DofVector]:
        """Loop over all elements, and evalute them.
        Is is called by solveStep() in each iteration.

        Parameters
        ----------
        elements
            The list of finite elements.
        U_np
            The current solution vector.
        dU
            The current solution increment vector.
        P
            The reaction vector.
        K
            The system matrix.
        F
            The vector of accumulated fluxes for convergence checks.
        increment
            The increment.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix,DofVector]
            - The modified reaction vector.
            - The modified system matrix.
            - The modified accumulated flux vector.
        """

        tic = getCurrentTime()

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        time = np.array([stepTime, totalTime])

        for el in elements.values():

            Ke = K[el]
            Pe = np.zeros(el.nDof)

            el.computeYourself(Ke, Pe, U_np[el], dU[el], time, dT)

            P[el] += Pe
            F[el] += abs(Pe)

        toc = getCurrentTime()
        self.computationTimes["elements"] += toc - tic

        return P, K, F

    def assembleConstraints(
        self,
        constraints: list[ConstraintBase],
        U_np: DofVector,
        dU: DofVector,
        PExt: DofVector,
        K: VIJSystemMatrix,
        increment: tuple,
    ) -> tuple[DofVector, VIJSystemMatrix]:
        """Loop over all elements, and evaluate them.
        Is is called by solveStep() in each iteration.

        Parameters
        ----------
        constraints
            The list of constraints.
        U_np
            The current solution vector.
        dU
            The current solution increment vector.
        PExt
            The external load vector.
        K
            The system matrix.
        increment
            The increment.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix,DofVector]
            - The modified external load vector.
            - The modified system matrix.
        """

        tic = getCurrentTime()

        for constraint in constraints.values():

            Ke = K[constraint]
            Pe = np.zeros(constraint.nDof)

            constraint.applyConstraint(U_np[constraint], dU[constraint], Pe, Ke, increment)

            # instead of PExt[constraint] += Pe, np.add.at allows for repeated indices
            np.add.at(PExt, PExt.entitiesInDofVector[constraint], Pe)

        toc = getCurrentTime()
        self.computationTimes["constraints"] += toc - tic

        return PExt, K

    def assembleLoads(
        self,
        nodeForces: list[StepActionBase],
        distributedLoads: list[StepActionBase],
        bodyForces: list[StepActionBase],
        U_np: DofVector,
        PExt: DofVector,
        K: VIJSystemMatrix,
        increment: tuple,
    ) -> tuple[DofVector, VIJSystemMatrix]:
        """Assemble all loads into a right hand side vector.

        Parameters
        ----------
        nodeForces
            The list of concentrated (nodal) loads.
        distributedLoads
            The list of distributed (surface) loads.
        bodyForces
            The list of body (volumetric) loads.
        U_np
            The current solution vector.
        PExt
            The external load vector.
        K
            The system matrix.
        increment
            The increment.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix]
            - The modified external load vector.
            - The modified system matrix.
        """
        # cloads
        for cLoad in nodeForces:
            PExt = cLoad.applyOnP(PExt, increment)
        # dloads
        PExt, K = self.computeDistributedLoads(distributedLoads, U_np, PExt, K, increment)
        PExt, K = self.computeBodyForces(bodyForces, U_np, PExt, K, increment)

        return PExt, K

    def extrapolateLastIncrement(
        self, extrapolation: str, increment: tuple, dU: DofVector, dirichlets: list, lastIncrementSize: float
    ) -> tuple[DofVector, bool]:
        """Depending on the current setting, extrapolate the solution of the last increment.

        Parameters
        ----------
        extrapolation
            The type of extrapolation.
        increment
            The current time increment.
        dU
            The last solution increment.
        dirichlets
            The list of active dirichlet boundary conditions.
        lastIncrementSize
            The size of the last increment.

        Returns
        -------
        tuple[DofVector,bool]
            - The extrapolated solution increment.
            - True if an extrapolation was performed.
        """

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment

        if extrapolation == "linear" and lastIncrementSize:
            dU *= incrementSize / lastIncrementSize
            dU = self.applyDirichlet(increment, dU, dirichlets)
            isExtrapolatedIncrement = True
        else:
            isExtrapolatedIncrement = False
            dU[:] = 0.0

        return dU, isExtrapolatedIncrement

    def checkDivergingSolution(self, incrementResidualHistory: dict, maxGrowingIter: int) -> bool:
        """Check if the iterative solution scheme is diverging.

        Parameters
        ----------
        incrementResidualHistory
            The dictionary containing the residual history of all fields.
        maxGrowingIter
            The maximum allows number of growths of a residual during the iterative solution scheme.

        Returns
        -------
        bool
            True if solution is diverging.
        """
        for previousFluxResidual, nGrew in incrementResidualHistory.values():
            if nGrew > maxGrowingIter:
                return True
        return False

    def printResidualOutlierNodes(self, residualOutliers: dict):
        """Print which nodes have the largest residuals.

        Parameters
        ----------
        residualOutliers
            The dictionary containing the outlier nodes for every field.
        """
        self.journal.message(
            "Residual outliers:",
            self.identification,
            level=1,
        )
        for field, node in residualOutliers.items():
            self.journal.message(
                "|{:20}|node {:10}|".format(field, node.label),
                self.identification,
                level=2,
            )

    def applyStepActionsAtStepStart(self, U: DofVector, P: DofVector, stepActions: dict[str, StepActionBase]):
        """Called when all step actions should be appliet at the start a step.

        Parameters
        ----------
        U
            The solution vector.
        P
            The reaction vector.
        stepActions
            The dictionary of active step actions.
        """

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtStepStart(U, P)

    def applyStepActionsAtStepEnd(self, U: DofVector, P: DofVector, stepActions: dict[str, StepActionBase]):
        """Called when all step actions should finish a step.

        Parameters
        ----------
        U
            The solution vector.
        P
            The reaction vector.
        stepActions
            The dictionary of active step actions.
        """

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtStepEnd(U, P)

    def applyStepActionsAtIncrementStart(
        self, U_n: DofVector, P: DofVector, increment: tuple, stepActions: dict[str, StepActionBase]
    ):
        """Called when all step actions should be applied at the start of a step.

        Parameters
        ----------
        U_n
            The converged solution vector at the strt of an increment.
        P
            The reaction vector.
        increment
            The time increment.
        stepActions
            The dictionary of active step actions.
        """

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtIncrementStart(U_n, P, increment)
