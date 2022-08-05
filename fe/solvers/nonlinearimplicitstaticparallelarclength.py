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


class NISTPArcLength(NISTParallel):
    identification = "NISTPArcLength"

    def __init__(self, jobInfo, journal):

        self.Lambda = 0.0

        return super().__init__(jobInfo, journal)

    def solveStep(self, step, time, stepActions, modelInfo, U, P, fieldOutputController, outputmanagers):

        self.arcLengthController = None
        self.checkConditionalStop = lambda: False
        # self.

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
            # if stopCondition != "False":
            self.checkConditionalStop = createModelAccessibleFunction(
                stopCondition, modelInfo=modelInfo, fieldOutputs=fieldOutputController.fieldOutputs
            )
        except KeyError:
            pass

        self.dLambda = 0.0

        # if "stopCondition" in options:
        # else:

        return super().solveStep(step, time, stepActions, modelInfo, U, P, fieldOutputController, outputmanagers)

    def solveIncrement(
        self,
        U,
        dU,
        P,
        K,
        elements,
        activeStepActions,
        constraints,
        increment,
        lastIncrementSize,
        extrapolation,
        maxIter,
        maxGrowingIter,
    ):
        """Arc length method to solve for an increment
        Implementation based on the proposed approach by"""

        if self.arcLengthController == None:
            return super().solveIncrement(
                U,
                dU,
                P,
                K,
                elements,
                activeStepActions,
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

        dirichlets = activeStepActions["dirichlets"]
        concentratedLoads = activeStepActions["concentratedLoads"]
        distributedLoads = activeStepActions["distributedLoads"]
        bodyForces = activeStepActions["bodyForces"]

        R_ = np.tile(self.theDofManager.constructDofVector(), (2, 1)).T  # 2 RHSs
        R_0 = R_[:, 0]
        R_f = R_[:, 1]
        F = self.theDofManager.constructDofVector()  # accumulated Flux vector

        P_0 = self.theDofManager.constructDofVector()
        P_f = self.theDofManager.constructDofVector()
        K_f = self.theDofManager.constructVIJSystemMatrix()
        K_0 = self.theDofManager.constructVIJSystemMatrix()
        Un1 = self.theDofManager.constructDofVector()
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
            for geostatic in activeStepActions["geostatics"]:
                geostatic.apply()

            Un1[:] = U
            Un1 += dU

            P[:] = K[:] = F[:] = P_0[:] = P_f[:] = K_f[:] = K_0[:] = 0.0

            P, K, F = self.computeElements(elements, Un1, dU, P, K, F, increment)
            P, K = self.assembleConstraints(constraints, Un1, dU, P, K, increment)

            P_0, K_0 = self.assembleLoads(
                concentratedLoads, distributedLoads, bodyForces, Un1, P_0, K_0, zeroIncrement
            )  # compute 'dead' deadloads, like gravity
            P_f, K_f = self.assembleLoads(
                concentratedLoads, distributedLoads, bodyForces, Un1, P_f, K_f, referenceIncrement
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
        self.arcLengthController.finishIncrement(U, dU, dLambda)

        self.journal.message(
            "Current load parameter: lambda={:6.2e}, last increment: dLambda={:6.2e}".format(self.Lambda, self.dLambda),
            self.identification,
            level=1,
        )

        return Un1, dU, P, iterationCounter, incrementResidualHistory

    def extrapolateLastIncrement(self, extrapolation, increment, dU, dirichlets, lastIncrementSize, dLambda=None):
        """Modified extrapolation to account for a load multiplier predictor"""

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

    def finishStepActions(self, U, P, stepActions):
        """Modified finish to communicate the 'correct' magnitude of external loads,
        i.e., the load multiplier,
        to the stepactions (nodeforces, distributedloads, ... )"""

        if self.arcLengthController != None:
            for stepAction in stepActions["nodeforces"].values():
                stepAction.finishStep(U, P, stepMagnitude=self.Lambda)
            for stepAction in stepActions["distributedload"].values():
                stepAction.finishStep(U, P, stepMagnitude=self.Lambda)

        return super().finishStepActions(U, P, stepActions)
