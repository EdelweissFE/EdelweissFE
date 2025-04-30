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

import json
from time import time as getCurrentTime

import numpy as np
from numpy import ndarray
from scipy.sparse import csr_matrix

from edelweissfe.config.linsolve import getDefaultLinSolver, getLinSolverByName
from edelweissfe.config.timing import createTimingDict
from edelweissfe.constraints.base.constraintbase import ConstraintBase
from edelweissfe.models.femodel import FEModel
from edelweissfe.numerics.csrgenerator import CSRGenerator
from edelweissfe.numerics.dofmanager import DofManager, DofVector, VIJSystemMatrix
from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from edelweissfe.stepactions.base.stepactionbase import StepActionBase
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.utils.exceptions import (
    ConditionalStop,
    CutbackRequest,
    DivergingSolution,
    ReachedMaxIncrements,
    ReachedMaxIterations,
    ReachedMinIncrementSize,
    StepFailed,
)
from edelweissfe.utils.fieldoutput import FieldOutputController


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

    NISTOptions = {
        "defaultMaxIter": 10,
        "defaultCriticalIter": 5,
        "defaultMaxGrowingIter": 10,
        "extrapolation": "linear",
        "linsolver": "pardiso",
        "linsolverConfigFile": "",
    }

    def __init__(self, jobInfo, journal, **kwargs):
        self.journal = journal

        self.fieldCorrectionTolerances = jobInfo["fieldCorrectionTolerance"]
        self.fluxResidualTolerances = jobInfo["fluxResidualTolerance"]
        self.fluxResidualTolerancesAlt = jobInfo["fluxResidualToleranceAlternative"]

        self.options = self.NISTOptions.copy()
        self._updateOptions(kwargs, journal)

    def _updateOptions(self, updatedOptions: dict, journal):
        """Update options of the solver using a string dict

        Parameters
        ----------
        updatedOptions
            The options dictionary.
        journal
            The journal module.
        """

        for k, v in updatedOptions.items():
            if k in self.NISTOptions:
                journal.message("Updating option {:}={:}".format(k, v), self.identification)
                self.options[k] = type(self.NISTOptions[k])(updatedOptions[k])
            else:
                raise AttributeError("Invalid option {:} for {:}".format(k, self.identification))

    def solveStep(
        self,
        step,
        model: FEModel,
        fieldOutputController: FieldOutputController,
        outputmanagers: dict[str, OutputManagerBase],
    ) -> tuple[bool, FEModel]:
        """Public interface to solve for a step.

        Parameters
        ----------
        stepNumber
            The step number.
        step
            The dictionary containing the step definition.
        stepActions
            The dictionary containing all step actions.
        model
            The  model tree.
        fieldOutputController
            The field output controller.
        """

        self.journal.message("Creating monolithic equation system", self.identification, 0)
        self.theDofManager = DofManager(
            model.nodeFields.values(),
            model.scalarVariables.values(),
            model.elements.values(),
            model.constraints.values(),
            model.nodeSets.values(),
        )
        self.journal.message(
            "total size of eq. system: {:}".format(self.theDofManager.nDof),
            self.identification,
            0,
        )

        self.journal.printSeperationLine()

        presentVariableNames = list(self.theDofManager.idcsOfFieldsInDofVector.keys())

        if self.theDofManager.idcsOfScalarVariablesInDofVector:
            presentVariableNames += [
                "scalar variables",
            ]

        nVariables = len(presentVariableNames)
        self.iterationHeader = ("{:^25}" * nVariables).format(*presentVariableNames)
        self.iterationHeader2 = (" {:<10}  {:<10}  ").format("||R||∞", "||ddU||∞") * nVariables
        self.iterationMessageTemplate = "{:11.2e}{:1}{:11.2e}{:1} "

        U = self.theDofManager.constructDofVector()
        K = self.theDofManager.constructVIJSystemMatrix()

        self.csrGenerator = CSRGenerator(K)

        self.computationTimes = createTimingDict()

        try:
            self._updateOptions(step.actions["options"]["NISTSolver"].options, self.journal)
        except KeyError:
            pass

        extrapolation = self.options["extrapolation"]
        linsolverOptions = self.options["linsolverConfigFile"]
        linsolverOptionDict = json.load(open(linsolverOptions, "r")) if linsolverOptions else ""
        self.linSolver = (
            getLinSolverByName(self.options["linsolver"], linsolverOptionDict)
            if "linsolver" in self.options
            else getDefaultLinSolver()
        )

        maxIter = step.maxIter
        criticalIter = step.criticalIter
        maxGrowingIter = step.maxGrowIter

        # nodes = model.nodes
        # elements = model.elements
        # constraints = model.constraints

        U = self.theDofManager.constructDofVector()
        P = self.theDofManager.constructDofVector()
        dU = self.theDofManager.constructDofVector()

        for fieldName, field in model.nodeFields.items():
            U = self.theDofManager.writeNodeFieldToDofVector(U, field, "U")

        for variable in model.scalarVariables.values():
            U[self.theDofManager.idcsOfScalarVariablesInDofVector[variable]] = variable.value

        prevTimeStep = None

        self.applyStepActionsAtStepStart(model, step.actions)

        try:
            for timeStep in step.getTimeStep():
                statusInfoDict = {
                    "step": step.number,
                    "inc": timeStep.number,
                    "iters": None,
                    "converged": False,
                    "time inc": timeStep.timeIncrement,
                    "time end": timeStep.totalTime,
                    "notes": "",
                }

                self.journal.printSeperationLine()
                self.journal.message(
                    "increment {:}: {:8f}, {:8f}; time {:10f} to {:10f}".format(
                        timeStep.number,
                        timeStep.stepProgressIncrement,
                        timeStep.stepProgress,
                        timeStep.totalTime - timeStep.timeIncrement,
                        timeStep.totalTime,
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
                        step.actions,
                        model,
                        timeStep,
                        prevTimeStep,
                        extrapolation,
                        maxIter,
                        maxGrowingIter,
                    )

                except CutbackRequest as e:
                    self.journal.message(str(e), self.identification, 1)
                    step.discardAndChangeIncrement(max(e.cutbackSize, 0.25))
                    prevTimeStep = None

                    statusInfoDict["iters"] = np.inf
                    statusInfoDict["notes"] = str(e)

                    tic = getCurrentTime()
                    for man in outputmanagers:
                        man.finalizeFailedIncrement(
                            statusInfoDict=statusInfoDict,
                            currentComputingTimes=self.computationTimes,
                        )
                    toc = getCurrentTime()
                    self.computationTimes["output"] += toc - tic

                except (ReachedMaxIterations, DivergingSolution) as e:
                    self.journal.message(str(e), self.identification, 1)
                    step.discardAndChangeIncrement(0.25)
                    prevTimeStep = None

                    statusInfoDict["iters"] = np.inf
                    statusInfoDict["notes"] = str(e)

                    tic = getCurrentTime()
                    for man in outputmanagers:
                        man.finalizeFailedIncrement(
                            statusInfoDict=statusInfoDict,
                            currentComputingTimes=self.computationTimes,
                        )
                    toc = getCurrentTime()
                    self.computationTimes["output"] += toc - tic

                else:
                    prevTimeStep = timeStep
                    if iterationCounter >= criticalIter:
                        step.preventIncrementIncrease()

                    # write results to nodes:
                    for fieldName, field in model.nodeFields.items():
                        self.theDofManager.writeDofVectorToNodeField(U, field, "U")
                        self.theDofManager.writeDofVectorToNodeField(P, field, "P")
                        self.theDofManager.writeDofVectorToNodeField(dU, field, "dU")

                    for variable in model.scalarVariables.values():
                        variable.value = U[self.theDofManager.idcsOfScalarVariablesInDofVector[variable]]

                    model.advanceToTime(timeStep.totalTime)

                    self.journal.message(
                        "Converged in {:} iteration(s)".format(iterationCounter),
                        self.identification,
                        1,
                    )

                    statusInfoDict["iters"] = iterationCounter
                    statusInfoDict["converged"] = True

                    fieldOutputController.finalizeIncrement()
                    tic = getCurrentTime()
                    for man in outputmanagers:
                        man.finalizeIncrement(
                            currentComputingTimes=self.computationTimes,
                            statusInfoDict=statusInfoDict,
                        )
                    toc = getCurrentTime()
                    self.computationTimes["output"] += toc - tic

        except (ReachedMaxIncrements, ReachedMinIncrementSize):
            self.journal.errorMessage("Incrementation failed", self.identification)
            raise StepFailed()

        except ConditionalStop:
            self.journal.message("Conditional Stop", self.identification)
            self.applyStepActionsAtStepEnd(model, step.actions)

        else:
            self.applyStepActionsAtStepEnd(model, step.actions)

        finally:
            self.journal.printTable(
                [("Time in {:}".format(k), " {:10.4e}s".format(v)) for k, v in self.computationTimes.items()],
                self.identification,
            )

    def solveIncrement(
        self,
        U_n: DofVector,
        dU: DofVector,
        P: DofVector,
        K: VIJSystemMatrix,
        stepActions: list,
        model: FEModel,
        timeStep: TimeStep,
        prevTimeStep: TimeStep,
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
        model
            The model tree.
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

        # incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = timeStep
        iterationCounter = 0
        incrementResidualHistory = dict.fromkeys(self.theDofManager.idcsOfFieldsInDofVector, (0.0, 0))

        elements = model.elements
        constraints = model.constraints

        R = self.theDofManager.constructDofVector()
        F = self.theDofManager.constructDofVector()
        PExt = self.theDofManager.constructDofVector()
        U_np = self.theDofManager.constructDofVector()
        ddU = None

        dirichlets = stepActions["dirichlet"].values()
        nodeforces = stepActions["nodeforces"].values()
        distributedLoads = stepActions["distributedload"].values()
        bodyForces = stepActions["bodyforce"].values()

        self.applyStepActionsAtIncrementStart(model, timeStep, stepActions)

        dU, isExtrapolatedIncrement = self.extrapolateLastIncrement(
            extrapolation, timeStep, dU, dirichlets, prevTimeStep, model
        )

        while True:
            for geostatic in stepActions["geostatic"].values():
                geostatic.applyAtIterationStart()

            U_np[:] = U_n
            U_np += dU

            P[:] = K[:] = F[:] = PExt[:] = 0.0

            P, K, F = self.computeElements(elements, U_np, dU, P, K, F, timeStep)
            PExt, K = self.assembleLoads(nodeforces, distributedLoads, bodyForces, U_np, PExt, K, timeStep)
            PExt, K = self.assembleConstraints(constraints, U_np, dU, PExt, K, timeStep)

            R[:] = P
            R += PExt

            if iterationCounter == 0 and not isExtrapolatedIncrement and dirichlets:
                # first iteration? apply dirichlet bcs and unconditionally solve
                R = self.applyDirichlet(timeStep, R, dirichlets)
            else:
                # iteration cycle 1 or higher, time to check the convergence
                for dirichlet in dirichlets:
                    R[self.findDirichletIndices(dirichlet)] = 0.0

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
        timeStep: TimeStep,
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
        timeStep
            The current time step.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix]
            The augmented load vector and system matrix.
        """

        tic = getCurrentTime()

        time = np.array([timeStep.stepTime, timeStep.totalTime])
        dT = timeStep.timeIncrement

        for dLoad in distributedLoads:
            load = dLoad.getCurrentLoad(timeStep)
            for faceID, elementSet in dLoad.surface.items():
                for el in elementSet:
                    Ke = K[el]
                    Pe = np.zeros(el.nDof)

                    el.computeDistributedLoad(dLoad.loadType, Pe, Ke, faceID, load, U_np[el], time, dT)

                    PExt[el] += Pe

        toc = getCurrentTime()
        self.computationTimes["distributed loads"] += toc - tic
        return PExt, K

    def computeBodyForces(
        self,
        bodyForces: list[StepActionBase],
        U_np: DofVector,
        PExt: DofVector,
        K: VIJSystemMatrix,
        timeStep: TimeStep,
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

        time = np.array([timeStep.stepTime, timeStep.totalTime])
        dT = timeStep.timeIncrement

        for bForce in bodyForces:
            force = bForce.getCurrentLoad(timeStep)
            for el in bForce.elementSet:
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

        if dirichlets:
            tic = getCurrentTime()
            for dirichlet in dirichlets:
                for row in self.findDirichletIndices(dirichlet):  # dirichlet.indices:
                    K.data[K.indptr[row] : K.indptr[row + 1]] = 0.0

            # K[row, row] = 1.0 @ once, faster than within the loop above:
            diag = K.diagonal()
            diag[np.concatenate([self.findDirichletIndices(d) for d in dirichlets])] = 1.0
            K.setdiag(diag)

            K.eliminate_zeros()

            toc = getCurrentTime()
            self.computationTimes["dirichlet K"] += toc - tic

        return K

    def applyDirichlet(self, timeStep: TimeStep, R: DofVector, dirichlets: list[StepActionBase]):
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
            delta = dirichlet.getDelta(timeStep)
            R[self.findDirichletIndices(dirichlet)] = delta.flatten()

        toc = getCurrentTime()
        self.computationTimes["dirichlet R"] += toc - tic

        return R

    def checkConvergence(
        self,
        R: DofVector,
        ddU: DofVector,
        F: DofVector,
        iterationCounter: int,
        residualHistory: dict,
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

        for field, fieldIndices in self.theDofManager.idcsOfFieldsInDofVector.items():
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
            convergedAtAll = convergedAtAll and convergedCorrection and convergedFlux

        if self.theDofManager.idcsOfScalarVariablesInDofVector:
            residualScalarVariables = max(np.abs(R[list(self.theDofManager.idcsOfScalarVariablesInDofVector.values())]))
            correction = (
                np.linalg.norm(
                    ddU[list(self.theDofManager.idcsOfScalarVariablesInDofVector.values())],
                    np.inf,
                )
                if ddU is not None
                else 0.0
            )

            convergedCorrection = correction < self.fieldCorrectionTolerances["scalar variables"]
            convergedFlux = residualScalarVariables <= fluxResidualTolerances["scalar variables"]

            iterationMessage += self.iterationMessageTemplate.format(
                residualScalarVariables,
                "✓" if convergedFlux else " ",
                correction,
                "✓" if convergedCorrection else " ",
            )

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

    def computeSpatialAveragedFluxes(self, F: DofVector) -> float:
        """Compute the spatial averaged flux for every field
        Is usually called by checkConvergence().

        Parameters
        ----------
        F
            The accumulated flux vector.

        Returns
        -------
        dict[str,float]
            A dictioary containg the spatial average fluxes for every field.
        """
        spatialAveragedFluxes = dict.fromkeys(self.theDofManager.idcsOfFieldsInDofVector, 0.0)
        for field, nDof in self.theDofManager.nAccumulatedNodalFluxesFieldwise.items():
            spatialAveragedFluxes[field] = max(
                1e-10,
                np.linalg.norm(F[self.theDofManager.idcsOfFieldsInDofVector[field]], 1) / nDof,
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
        timeStep: TimeStep,
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
        timeStep
            The time step.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix,DofVector]
            - The modified reaction vector.
            - The modified system matrix.
            - The modified accumulated flux vector.
        """

        tic = getCurrentTime()

        time = np.array([timeStep.stepTime, timeStep.totalTime])
        dT = timeStep.timeIncrement

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
        timeStep: TimeStep,
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
        dT
            The time increment.
        time
            The step and total time.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix,DofVector]
            - The modified external load vector.
            - The modified system matrix.
        """

        tic = getCurrentTime()

        for constraint in constraints.values():
            Kc = K[constraint].reshape(constraint.nDof, constraint.nDof, order="F")
            Pc = np.zeros(constraint.nDof)

            constraint.applyConstraint(U_np[constraint], dU[constraint], Pc, Kc, timeStep)

            # instead of PExt[constraint] += Pe, np.add.at allows for repeated indices
            np.add.at(PExt, PExt.entitiesInDofVector[constraint], Pc)

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
        timeStep: TimeStep,
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
        timeStep
            The current time step.

        Returns
        -------
        tuple[DofVector,VIJSystemMatrix]
            - The augmented external load vector.
            - The augmented system matrix.
        """
        for cLoad in nodeForces:
            PExt[
                self.theDofManager.idcsOfFieldsOnNodeSetsInDofVector[cLoad.field][cLoad.nodeSet]
            ] += cLoad.getCurrentLoad(timeStep).flatten()
        PExt, K = self.computeDistributedLoads(distributedLoads, U_np, PExt, K, timeStep)
        PExt, K = self.computeBodyForces(bodyForces, U_np, PExt, K, timeStep)

        return PExt, K

    def extrapolateLastIncrement(
        self,
        extrapolation: str,
        timeStep: TimeStep,
        dU: DofVector,
        dirichlets: list,
        prevTimeStep: TimeStep,
        model,
    ) -> tuple[DofVector, bool]:
        """Depending on the current setting, extrapolate the solution of the last increment.

        Parameters
        ----------
        extrapolation
            The type of extrapolation.
        timeStep
            The current time step.
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

        if extrapolation == "linear" and prevTimeStep and prevTimeStep.timeIncrement:
            dU *= timeStep.stepProgressIncrement / prevTimeStep.stepProgressIncrement
            dU = self.applyDirichlet(timeStep, dU, dirichlets)
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

    def applyStepActionsAtStepStart(self, model: FEModel, stepActions: dict[str, StepActionBase]):
        """Called when all step actions should be appliet at the start a step.

        Parameters
        ----------
        model
            The model tree.
        stepActions
            The dictionary of active step actions.
        """

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtStepStart(model)

    def applyStepActionsAtStepEnd(self, model: FEModel, stepActions: dict[str, StepActionBase]):
        """Called when all step actions should finish a step.

        Parameters
        ----------
        model
            The model tree.
        stepActions
            The dictionary of active step actions.
        """

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtStepEnd(model)

    def applyStepActionsAtIncrementStart(
        self, model: FEModel, timeStep: TimeStep, stepActions: dict[str, StepActionBase]
    ):
        """Called when all step actions should be applied at the start of a step.

        Parameters
        ----------
        model
            The model tree.
        increment
            The time increment.
        stepActions
            The dictionary of active step actions.
        """

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtIncrementStart(model, timeStep)

    def findDirichletIndices(self, dirichlet):
        nSet = dirichlet.nSet
        field = dirichlet.field
        components = dirichlet.components

        fieldIndices = self.theDofManager.idcsOfFieldsOnNodeSetsInDofVector[field][nSet]

        return fieldIndices.reshape((len(nSet), -1))[:, components].flatten()
