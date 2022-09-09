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
# Created on Thu Nov  2 13:43:01 2017

# @author: Matthias Neuner
"""
(Parallel) Arc Length Solver, based on the proposed approach in Jirásek/Bažant 2001.
Replaces the NewtonRaphson scheme of the NISTParallel Solver.
"""

from fe.solvers.nonlinearimplicitstaticparallel import NISTParallel

import numpy as np
from fe.utils.exceptions import ReachedMaxIterations, DivergingSolution, ConditionalStop
from fe.utils.math import createModelAccessibleFunction
from fe.utils.dofmanager import DofVector, VIJSystemMatrix
from fe.stepactions.base.stepactionbase import StepActionBase
from fe.utils.fieldoutput import FieldOutputController
from fe.outputmanagers.base.outputmanagerbase import OutputManagerBase


class NISTPArcLength(NISTParallel):
    identification = "NISTPArcLength"

    def __init__(self, jobInfo, journal):

        self.Lambda = 0.0

        return super().__init__(jobInfo, journal)

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

        self.arcLengthController = None
        self.checkConditionalStop = lambda: False

        try:
            arcLengthControllerType = stepActions["options"]["NISTArcLength"]["arcLengthController"]
            if arcLengthControllerType == "off":
                self.arcLengthController = None
            else:
                self.arcLengthController = stepActions[arcLengthControllerType][arcLengthControllerType]
        except KeyError:
            pass

        try:
            stopCondition = stepActions["options"]["NISTArcLength"]["stopCondition"]
            self.checkConditionalStop = createModelAccessibleFunction(
                stopCondition, modelInfo=modelInfo, fieldOutputs=fieldOutputController.fieldOutputs
            )
        except KeyError:
            pass

        self.dLambda = 0.0

        return super().solveStep(
            stepNumber, step, time, stepActions, modelInfo, U, P, fieldOutputController, outputmanagers
        )

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
        """Arc length method scheme to solve for an increment.

        Parameters
        ----------
        U_n
            The old solution vector.
        dU
            The old solution increment.
        P
            The old reaction vector.
        K
            The system matrix to be used.
        elements
            The dictionary containing all elements.
        activeStepActions
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

        if self.arcLengthController == None:
            return super().solveIncrement(
                U_n,
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

        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment

        if incNumber > 1 and self.checkConditionalStop():
            raise ConditionalStop

        iterationCounter = 0
        incrementResidualHistory = dict.fromkeys(self.theDofManager.indicesOfFieldsInDofVector, (0.0, 0))

        dirichlets = stepActions["dirichlet"].values()
        nodeForces = stepActions["nodeforces"].values()
        distributedLoads = stepActions["distributedload"].values()
        bodyForces = stepActions["bodyforce"].values()

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtIncrementStart(U_n, P, increment)

        R_ = np.tile(self.theDofManager.constructDofVector(), (2, 1)).T  # 2 RHSs
        R_0 = R_[:, 0]
        R_f = R_[:, 1]
        F = self.theDofManager.constructDofVector()  # accumulated Flux vector

        P_0 = self.theDofManager.constructDofVector()
        P_f = self.theDofManager.constructDofVector()
        K_f = self.theDofManager.constructVIJSystemMatrix()
        K_0 = self.theDofManager.constructVIJSystemMatrix()
        U_np = self.theDofManager.constructDofVector()
        ddU = None

        Lambda = self.Lambda
        dLambda = self.dLambda
        ddLambda = 0.0

        dU, isExtrapolatedIncrement, dLambda = self.extrapolateLastIncrement(
            extrapolation, increment, dU, dirichlets, lastIncrementSize, dLambda
        )

        referenceIncrement = incNumber, 1.0, 1.0, 0.0, 0.0, 0.0
        zeroIncrement = incNumber, 0.0, 0.0, 0.0, 0.0, 0.0

        while True:
            for geostatic in stepActions["geostatics"]:
                geostatic.apply()

            U_np[:] = U_n
            U_np += dU

            P[:] = K[:] = F[:] = P_0[:] = P_f[:] = K_f[:] = K_0[:] = 0.0

            P, K, F = self.computeElements(elements, U_np, dU, P, K, F, increment)
            P, K = self.assembleConstraints(constraints, U_np, dU, P, K, increment)

            P_0, K_0 = self.assembleLoads(
                nodeForces, distributedLoads, bodyForces, U_np, P_0, K_0, zeroIncrement
            )  # compute 'dead' deadloads, like gravity
            P_f, K_f = self.assembleLoads(
                nodeForces, distributedLoads, bodyForces, U_np, P_f, K_f, referenceIncrement
            )  # compute 'dead' deadloads, like gravity

            P_f -= P_0  # and subtract the dead part, since we are only interested in the homogeneous linear part
            K_f -= K_0

            # Dead and Reference load ..
            R_0[:] = P_0 + (Lambda + dLambda) * P_f + P
            R_f[:] = P_f

            # add stiffness contribution
            K[:] += K_0
            K[:] += (Lambda + dLambda) * K_f

            # Dirichlets ..
            if isExtrapolatedIncrement and iterationCounter == 0:
                R_0 = self.applyDirichlet(zeroIncrement, R_0, dirichlets)
            else:
                modifiedIncrement = incNumber, dLambda, Lambda + dLambda, 0.0, 0.0, 0.0
                R_0 = self.applyDirichlet(modifiedIncrement, R_0, dirichlets)

            R_f = self.applyDirichlet(referenceIncrement, R_f, dirichlets)

            if iterationCounter > 0 or isExtrapolatedIncrement:
                converged, nodesWithLargestResidual = self.checkConvergence(
                    R_0, ddU, F, iterationCounter, incrementResidualHistory
                )

                if converged:
                    break

                if self.checkDivergingSolution(incrementResidualHistory, maxGrowingIter):
                    self.printResidualOutlierNodes(nodesWithLargestResidual)
                    raise DivergingSolution("Residual grew {:} times, cutting back".format(maxGrowingIter))

                if iterationCounter == maxIter:
                    self.printResidualOutlierNodes(nodesWithLargestResidual)
                    raise ReachedMaxIterations("Reached max. iterations in current increment, cutting back")

            # if iterationCounter == maxIter:
            #     raise ReachedMaxIterations("Reached max. iterations in current increment, cutting back")

            K_ = self.assembleStiffnessCSR(K)
            K_ = self.applyDirichletK(K_, dirichlets)

            # solve 2 eq. systems at once:
            ddU_ = self.linearSolve(K_, R_)
            # q_0 = K⁻¹ * (  Pext_0  + dLambda * Pext_Ref - PInt  )
            # q_f = K⁻¹ * (  Pext_Ref  )
            ddU_0, ddU_f = ddU_[:, 0], ddU_[:, 1]

            # compute the increment of the load parameter. Method depends on the employed arc length controller
            ddLambda = self.arcLengthController.computeDDLambda(dU, ddU_0, ddU_f, increment)

            # assemble total solution
            ddU = ddU_0 + ddLambda * ddU_f

            dU += ddU
            dLambda += ddLambda

            iterationCounter += 1

        self.Lambda += dLambda
        self.dLambda = dLambda
        self.arcLengthController.finishIncrement(U_n, dU, dLambda)

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtIncrementStart(U_n, P, increment)

        self.journal.message(
            "Current load parameter: lambda={:6.2e}, last increment: dLambda={:6.2e}".format(self.Lambda, self.dLambda),
            self.identification,
            level=1,
        )

        return U_np, dU, P, iterationCounter, incrementResidualHistory

    def extrapolateLastIncrement(
        self,
        extrapolation: str,
        increment: tuple,
        dU: DofVector,
        dirichlets: list,
        lastIncrementSize: float,
        dLambda: float = None,
    ) -> tuple[DofVector, bool]:
        """Depending on the current setting, extrapolate the solution of the last increment.
        Also accounts for the arc length loading parameter.

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

        if dLambda == None:
            # arclength control is inactive
            return super().extrapolateLastIncrement(extrapolation, increment, dU, dirichlets, lastIncrementSize)

        elif extrapolation == "linear" and lastIncrementSize:
            dLambda = dLambda * (incrementSize / lastIncrementSize)
            dU, isExtrapolatedIncrement = super().extrapolateLastIncrement(
                extrapolation, increment, dU, {}, lastIncrementSize
            )
        else:
            dLambda = 0.0
            dU[:] = 0.0
            isExtrapolatedIncrement = False

        return dU, isExtrapolatedIncrement, dLambda

    def applyStepActionsAtStepEnd(self, U: DofVector, P: DofVector, stepActions: dict[str, StepActionBase]):
        """Called when all step actions should finish a step.

        For the arc length controlled simulation,
        external loads at the end of a step depend on the arc length parameter.
        Hence, for subsequent steps, the new reference loads depend on this parameter.
        We communicate this to the step actions before the parent solver tells them that the
        step is finished.

        Parameters
        ----------
        U
            The solution vector.
        P
            The reaction vector.
        stepActions
            The list of active step actions.
        """

        if self.arcLengthController != None:
            for stepAction in stepActions["nodeforces"].values():
                stepAction.applyAtStepEnd(U, P, stepMagnitude=self.Lambda)
            for stepAction in stepActions["distributedload"].values():
                stepAction.applyAtStepEnd(U, P, stepMagnitude=self.Lambda)
            for stepAction in stepActions["bodeforce"].values():
                stepAction.applyAtStepEnd(U, P, stepMagnitude=self.Lambda)

        return super().applyStepActionsAtStepEnd(U, P, stepActions)
