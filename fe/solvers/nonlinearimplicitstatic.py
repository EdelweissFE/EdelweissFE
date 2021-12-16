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
"""
Created on Sun Jan  8 20:37:35 2017

@author: Matthias Neuner

Standard nonlinear, implicit static solver.

"""
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


class NIST:
    """This is the Nonlinear Implicit STatic -- solver.
    Designed to interface with Abaqus UELs
    Public methods are: __init__(), initializeUP() and solveStep(...).
    OutputManagers are updated at the end of each increment.
    """

    identification = "NISTSolver"

    defaultMaxInc = 1.0
    defaultMinInc = 1e-4
    defaultMaxNumInc = 1000
    defaultMaxIter = 10
    defaultCriticalIter = 5
    defaultMaxGrowingIter = 10

    def __init__(self, jobInfo, modelInfo, journal, fieldOutputController, outputmanagers):
        self.nodes = modelInfo["nodes"]
        self.elements = modelInfo["elements"]
        self.nodeSets = modelInfo["nodeSets"]
        self.elementSets = modelInfo["elementSets"]
        self.constraints = modelInfo["constraints"]

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
        self.fieldOutputController = fieldOutputController
        self.outputmanagers = outputmanagers

        self.systemMatrix = self.theDofManager.constructVIJSystemMatrix()

        self.csrGenerator = CSRGenerator(self.systemMatrix)

        self.residualHistories = dict.fromkeys(self.theDofManager.indicesOfFieldsInDofVector)

    def initialize(self):
        """Initialize the solver and return the 2 vectors for field (U) and flux (P)"""

        U = self.theDofManager.constructDofVector()
        P = self.theDofManager.constructDofVector()

        return U, P

    def solveStep(self, step, time, stepActions, stepOptions, U, P):
        """Public interface to solve for an ABAQUS like step
        returns: boolean Success, U vector, P vector, and the new current total time"""

        self.computationTimes = defaultdict(lambda: 0.0)

        K = self.systemMatrix

        stepLength = step.get("stepLength", 1.0)
        extrapolation = stepOptions["NISTSolver"].get("extrapolation", "linear")

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

        activeStepActions = self.collectActiveStepActions(stepActions)

        lastIncrementSize = False

        for setfield in activeStepActions["setfields"]:
            setfield.apply()
        for initmaterial in activeStepActions["initializematerial"]:
            initmaterial.apply()

        try:
            for increment in incGen.generateIncrement():
                incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment

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
                        activeStepActions,
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

                except (ReachedMaxIterations, DivergingSolution) as e:
                    self.journal.message(str(e), self.identification, 1)
                    incGen.discardAndChangeIncrement(0.25)
                    lastIncrementSize = False

                else:
                    lastIncrementSize = incrementSize
                    if iterationCounter >= criticalIter:
                        incGen.preventIncrementIncrease()

                    for el in self.elements.values():
                        el.acceptLastState()

                    self.journal.message(
                        "Converged in {:} iteration(s)".format(iterationCounter), self.identification, 1
                    )

                    self.fieldOutputController.finalizeIncrement(U, P, increment)
                    for man in self.outputmanagers:
                        man.finalizeIncrement(U, P, increment)

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

            self.finishStepActions(U, P, stepActions)

        else:
            success = True

            if activeStepActions["geostatics"]:
                self.journal.message("Geostatic Step -- displacements are resetted", self.identification)
                U = self.resetDisplacements(U)  # reset all displacements, if the present step is a geostatic step

            self.finishStepActions(U, P, stepActions)

        finally:
            finishedTime = time + stepProgress * stepLength
            self.journal.printTable(
                [("Time in {:}".format(k), " {:10.4f}s".format(v)) for k, v in self.computationTimes.items()],
                self.identification,
            )

        return success, U, P, finishedTime

    def solveIncrement(
        self, Un, dU, P, K, activeStepActions, increment, lastIncrementSize, extrapolation, maxIter, maxGrowingIter
    ):
        """Standard Newton-Raphson scheme to solve for an increment"""

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        iterationCounter = 0
        incrementResidualHistory = dict.fromkeys(self.theDofManager.indicesOfFieldsInDofVector, (0.0, 0))

        R = self.theDofManager.constructDofVector()
        F = self.theDofManager.constructDofVector()
        PExt = self.theDofManager.constructDofVector()
        Un1 = self.theDofManager.constructDofVector()
        ddU = None

        dirichlets = activeStepActions["dirichlets"]
        concentratedLoads = activeStepActions["concentratedLoads"]
        distributedLoads = activeStepActions["distributedLoads"]
        bodyForces = activeStepActions["bodyForces"]
        constraints = activeStepActions["constraints"]

        dU, isExtrapolatedIncrement = self.extrapolateLastIncrement(
            extrapolation, increment, dU, dirichlets, lastIncrementSize
        )

        while True:
            for geostatic in activeStepActions["geostatics"]:
                geostatic.apply()

            Un1[:] = Un
            Un1 += dU

            P[:] = K[:] = F[:] = PExt[:] = 0.0

            P, K, F = self.computeElements(Un1, dU, P, K, F, increment)
            PExt, K = self.assembleLoads(concentratedLoads, distributedLoads, bodyForces, Un1, PExt, K, increment)
            PExt, K = self.assembleConstraints(constraints, Un1, PExt, K, increment)

            R[:] = P
            R += PExt

            if iterationCounter == 0 and not isExtrapolatedIncrement and dirichlets:
                # first iteration? apply dirichlet bcs and unconditionally solve
                R = self.applyDirichlet(increment, R, dirichlets)
            else:
                # iteration cycle 1 or higher, time to check the convergence
                for dirichlet in dirichlets:
                    R[dirichlet.indices] = 0.0
                if self.checkConvergence(R, ddU, F, iterationCounter, incrementResidualHistory):
                    break

                if self.checkDivergingSolution(incrementResidualHistory, maxGrowingIter):
                    raise DivergingSolution("Residual grew {:} times, cutting back".format(maxGrowingIter))

            if iterationCounter == maxIter:
                raise ReachedMaxIterations("Reached max. iterations in current increment, cutting back")

            K_ = self.assembleStiffnessCSR(K)
            K_ = self.applyDirichletK(K_, dirichlets)

            ddU = self.linearSolve(K_, R)
            dU += ddU
            iterationCounter += 1

        return Un1, dU, P, iterationCounter, incrementResidualHistory

    def computeElements(self, Un1, dU, P, K, F, increment):
        """Loop over all elements, and evalute them.
        Note that ABAQUS style is employed: element(Un+1, dUn+1)
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration"""

        tic = getCurrentTime()

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        time = np.array([stepTime, totalTime])

        for el in self.elements.values():

            Ke = K[el]
            Pe = np.zeros(el.nDof)

            el.computeYourself(Ke, Pe, Un1[el], dU[el], time, dT)

            P[el] += Pe
            F[el] += abs(Pe)

        toc = getCurrentTime()
        self.computationTimes["elements"] += toc - tic

        return P, K, F

    def computeDistributedLoads(self, distributedLoads, Un1, PExt, K, increment):
        """Loop over all elements, and evalute them.
        Note that ABAQUS style is employed: element(Un+1, dUn+1)
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration"""

        tic = getCurrentTime()
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        time = np.array([stepTime, totalTime])
        for dLoad in distributedLoads:
            magnitude = dLoad.getCurrentMagnitude(increment)
            for faceID, elements in dLoad.surface.items():
                for el in elements:

                    Ke = K[el]
                    Pe = np.zeros(el.nDof)

                    el.computeDistributedLoad(dLoad.loadType, Pe, Ke, faceID, magnitude, Un1[el], time, dT)

                    PExt[el] += Pe

        toc = getCurrentTime()
        self.computationTimes["distributed loads"] += toc - tic
        return PExt, K

    def computeBodyForces(self, bodyForces, Un1, PExt, K, increment):
        """Loop over all elements, and evalute them.
        Note that ABAQUS style is employed: element(Un+1, dUn+1)
        instead of element(Un, dUn+1)
        -> is called by solveStep() in each iteration"""

        tic = getCurrentTime()
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        time = np.array([stepTime, totalTime])
        for bForce in bodyForces:
            force = bForce.getCurrentBodyForce(increment)
            for el in bForce.elements:

                Pe = np.zeros(el.nDof)
                Ke = K[el]

                el.computeBodyForce(Pe, Ke, force, Un1[el], time, dT)

                PExt[el] += Pe

        toc = getCurrentTime()
        self.computationTimes["body forces"] += toc - tic
        return PExt, K

    def applyDirichletK(self, K, dirichlets):
        """Apply the dirichlet bcs on the global stiffnes matrix
        -> is called by solveStep() before solving the global sys.
        http://stackoverflux.com/questions/12129948/scipy-sparse-set-row-to-zeros"""

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

    def applyDirichlet(self, increment, R, dirichlets):
        """Apply the dirichlet bcs on the Residual vector
        -> is called by solveStep() before solving the global sys."""

        tic = getCurrentTime()
        for dirichlet in dirichlets:
            R[dirichlet.indices] = dirichlet.getDelta(increment)

        toc = getCurrentTime()
        self.computationTimes["dirichlet R"] += toc - tic

        return R

    def checkConvergence(self, R, ddU, F, iterationCounter, residualHistory):
        """Check the convergency, indivudually for each field,
        similar to ABAQUS based on the current total flux residual and the field correction
        -> is called by solveStep() to decide wether to continue iterating or stop"""

        tic = getCurrentTime()

        iterationMessage = ""
        convergedAtAll = True

        spatialAveragedFluxes = self.computeSpatialAveragedFluxes(F)

        if iterationCounter < 15:  # standard tolerance set
            fluxResidualTolerances = self.fluxResidualTolerances
        else:  # alternative tolerance set
            fluxResidualTolerances = self.fluxResidualTolerancesAlt

        for field, fieldIndices in self.theDofManager.indicesOfFieldsInDofVector.items():
            fluxResidual = np.linalg.norm(R[fieldIndices], np.inf)
            fieldCorrection = np.linalg.norm(ddU[fieldIndices], np.inf) if ddU is not None else 0.0

            convergedCorrection = fieldCorrection < self.fieldCorrectionTolerances[field]
            convergedFlux = fluxResidual <= max(fluxResidualTolerances[field] * spatialAveragedFluxes[field], 1e-10)

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

            # converged if residual and fieldCorrection are smaller than tolerance
            convergedAtAll = convergedAtAll and convergedCorrection and convergedFlux

        self.journal.message(iterationMessage, self.identification)

        toc = getCurrentTime()
        self.computationTimes["convergence check"] += toc - tic

        return convergedAtAll

    def linearSolve(self, A, b):
        """Interface to the linear eq. system solver"""

        tic = getCurrentTime()
        ddU = self.linSolver(A, b)
        toc = getCurrentTime()
        self.computationTimes["linear solve"] += toc - tic
        return ddU

    def assembleStiffnessCSR(self, K):
        """Construct a CSR matrix from VIJ"""
        tic = getCurrentTime()
        KCsr = self.csrGenerator.updateCSR(K)
        toc = getCurrentTime()
        self.computationTimes["CSR generation"] += toc - tic
        return KCsr

    def resetDisplacements(self, U):
        """-> May be called at the end of a geostatic step, the zero all displacments after
        applying the geostatic stress"""

        self.journal.printSeperationLine()
        self.journal.message("Geostatic step, resetting displacements", self.identification, level=1)
        self.journal.printSeperationLine()
        U[self.theDofManager.indicesOfFieldsInDofVector["displacement"]] = 0.0
        return U

    def computeSpatialAveragedFluxes(self, F):
        """Compute the spatial averaged flux for every field
        -> Called by checkConvergence()"""
        spatialAveragedFluxes = dict.fromkeys(self.theDofManager.indicesOfFieldsInDofVector, 0.0)
        for field, nDof in self.theDofManager.nAccumulatedNodalFluxesFieldwise.items():
            spatialAveragedFluxes[field] = max(
                1e-10, np.linalg.norm(F[self.theDofManager.indicesOfFieldsInDofVector[field]], 1) / nDof
            )

        return spatialAveragedFluxes

    def assembleConstraints(self, constraints, Un1, PExt, K, increment):
        """Assemble the sub stiffness matrix for a linear constraint (independent of the solution)
        -> once per Increment"""

        tic = getCurrentTime()

        for constraint in constraints:

            Ke = K[constraint]
            Pe = np.zeros(constraint.nDof)

            constraint.applyConstraint(Un1[constraint], Pe, Ke, increment)

            PExt[constraint] += Pe

        toc = getCurrentTime()
        self.computationTimes["constraints"] += toc - tic

        return PExt, K

    def assembleLoads(self, concentratedLoads, distributedLoads, bodyForces, Un1, PExt, K, increment):
        """Assemble all loads into a right hand side vector"""
        # cloads
        for cLoad in concentratedLoads:
            PExt = cLoad.applyOnP(PExt, increment)
        # dloads
        PExt, K = self.computeDistributedLoads(distributedLoads, Un1, PExt, K, increment)
        PExt, K = self.computeBodyForces(bodyForces, Un1, PExt, K, increment)

        return PExt, K

    def extrapolateLastIncrement(self, extrapolation, increment, dU, dirichlets, lastIncrementSize):
        """if active, extrapolate the solution of the last increment.
        Also returns a boolean for information of an extrapolation was computed
        -> at the beginning of each increment"""
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment

        if extrapolation == "linear" and lastIncrementSize:
            dU *= incrementSize / lastIncrementSize
            dU = self.applyDirichlet(increment, dU, dirichlets)
            isExtrapolatedIncrement = True
        else:
            isExtrapolatedIncrement = False
            dU[:] = 0.0

        return dU, isExtrapolatedIncrement

    def checkDivergingSolution(self, incrementResidualHistory, maxGrowingIter):
        for previousFluxResidual, nGrew in incrementResidualHistory.values():
            if nGrew > maxGrowingIter:
                return True
        return False

    def collectActiveStepActions(self, stepActions):

        activeActions = dict()

        activeActions["dirichlets"] = list(stepActions["dirichlet"].values()) + list(
            stepActions["shotcreteshellmaster"].values()
        )
        activeActions["distributedLoads"] = stepActions["distributedload"].values()
        activeActions["bodyForces"] = stepActions["bodyforce"].values()
        activeActions["concentratedLoads"] = stepActions["nodeforces"].values()

        activeActions["geostatics"] = [g for g in stepActions["geostatic"].values() if g.active]
        activeActions["constraints"] = [c for c in self.constraints.values()]
        activeActions["setfields"] = [s for s in stepActions["setfield"].values() if s.active]
        activeActions["initializematerial"] = [s for s in stepActions["initializematerial"].values() if s.active]

        return activeActions

    def finishStepActions(self, U, P, stepActions):
        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.finishStep(U, P)
