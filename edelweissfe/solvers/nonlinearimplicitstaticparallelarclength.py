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

import numpy as np

from edelweissfe.models.femodel import FEModel
from edelweissfe.numerics.dofmanager import DofVector, VIJSystemMatrix
from edelweissfe.outputmanagers.base.outputmanagerbase import OutputManagerBase
from edelweissfe.solvers.nonlinearimplicitstaticparallel import NISTParallel
from edelweissfe.stepactions.base.stepactionbase import StepActionBase
from edelweissfe.timesteppers.timestep import TimeStep
from edelweissfe.utils.exceptions import (
    ConditionalStop,
    DivergingSolution,
    ReachedMaxIterations,
)
from edelweissfe.utils.fieldoutput import FieldOutputController
from edelweissfe.utils.math import createModelAccessibleFunction


class NISTPArcLength(NISTParallel):
    identification = "NISTPArcLength"

    def __init__(self, jobInfo, journal, **kwargs):
        self.Lambda = 0.0
        self.dLambda = 0.0
        self.arcLengthController = None

        return super().__init__(jobInfo, journal, **kwargs)

    def solveStep(
        self,
        step: dict,
        model: FEModel,
        fieldOutputController: FieldOutputController,
        outputmanagers: dict[str, OutputManagerBase],
    ):
        self.arcLengthController = None
        self.checkConditionalStop = lambda: False

        if "arc length parameter" in model.additionalParameters:
            self.Lambda = model.additionalParameters["arc length parameter"]

        arcLengthControllerOptions = step.actions["options"].get("NISTArcLength")
        if arcLengthControllerOptions:
            arcLengthController = arcLengthControllerOptions.get("arcLengthController")
            if arcLengthController:
                try:
                    self.arcLengthController = step.actions[arcLengthController][arcLengthController]
                    self.dLambda = 0.0
                except KeyError:
                    self.journal.errorMessage(
                        f'Arc length controller "{arcLengthController}" not found in current step or not configured correctly',
                        self.identification,
                    )
                    raise KeyError
            else:
                self.journal.message(
                    "No ArcLengthController specified in current step",
                    self.identification,
                )
                self.arcLengthController = None
                self.dLambda = None

            stopCondition = arcLengthControllerOptions.get("stopCondition")
            if stopCondition:
                # self.journal.message("St", self.identification)
                self.checkConditionalStop = createModelAccessibleFunction(
                    stopCondition,
                    model=model,
                    fieldOutputs=fieldOutputController.fieldOutputs,
                )
        else:
            self.journal.message("No ArcLengthController specified in current step", self.identification)
            self.arcLengthController = None
            self.dLambda = None

        return super().solveStep(step, model, fieldOutputController, outputmanagers)

    def solveIncrement(
        self,
        U_n: DofVector,
        dU: DofVector,
        P: DofVector,
        K: VIJSystemMatrix,
        stepActions: list,
        model,
        timeStep: TimeStep,
        prevTimeStep: TimeStep,
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
        timeStep
            The current time step.
        prevTimeStep
            The previous time step.
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

        if self.arcLengthController is None:
            return super().solveIncrement(
                U_n,
                dU,
                P,
                K,
                stepActions,
                model,
                timeStep,
                prevTimeStep,
                extrapolation,
                maxIter,
                maxGrowingIter,
            )

        if timeStep.number > 1 and self.checkConditionalStop():
            raise ConditionalStop

        iterationCounter = 0
        incrementResidualHistory = dict.fromkeys(self.theDofManager.idcsOfFieldsInDofVector, (0.0, 0))

        dirichlets = stepActions["dirichlet"].values()
        nodeForces = stepActions["nodeforces"].values()
        distributedLoads = stepActions["distributedload"].values()
        bodyForces = stepActions["bodyforce"].values()

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtIncrementStart(model, timeStep)

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
            extrapolation, timeStep, dU, dirichlets, prevTimeStep, model, dLambda
        )

        # referenceIncrement = incNumber, 1.0, 1.0, 0.0, 0.0, 0.0
        # zeroIncrement = incNumber, 0.0, 0.0, 0.0, 0.0, 0.0
        referenceTimeStep = TimeStep(timeStep.number, 1.0, 1.0, 0.0, 0.0, 0.0)
        zeroTimeStep = TimeStep(timeStep.number, 0.0, 0.0, 0.0, 0.0, 0.0)

        while True:
            for geostatic in stepActions["geostatics"]:
                geostatic.apply()

            U_np[:] = U_n
            U_np += dU

            P[:] = K[:] = F[:] = P_0[:] = P_f[:] = K_f[:] = K_0[:] = 0.0

            P, K, F = self.computeElements(model.elements, U_np, dU, P, K, F, timeStep)
            P, K = self.assembleConstraints(model.constraints, U_np, dU, P, K, timeStep)

            P_0, K_0 = self.assembleLoads(
                nodeForces, distributedLoads, bodyForces, U_np, P_0, K_0, zeroTimeStep
            )  # compute 'dead' deadloads, like gravity
            P_f, K_f = self.assembleLoads(
                nodeForces,
                distributedLoads,
                bodyForces,
                U_np,
                P_f,
                K_f,
                referenceTimeStep,
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
                R_0 = self.applyDirichlet(zeroTimeStep, R_0, dirichlets)
            else:
                modifiedTimeStep = TimeStep(timeStep.number, dLambda, Lambda + dLambda, 0.0, 0.0, 0.0)
                R_0 = self.applyDirichlet(modifiedTimeStep, R_0, dirichlets)

            R_f = self.applyDirichlet(referenceTimeStep, R_f, dirichlets)

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

            K_ = self.assembleStiffnessCSR(K)
            K_ = self.applyDirichletK(K_, dirichlets)

            # solve 2 eq. systems at once:
            ddU_ = self.linearSolve(K_, R_)
            # q_0 = K⁻¹ * (  Pext_0  + dLambda * Pext_Ref - PInt  )
            # q_f = K⁻¹ * (  Pext_Ref  )
            ddU_0, ddU_f = ddU_[:, 0], ddU_[:, 1]

            # compute the increment of the load parameter. Method depends on the employed arc length controller
            ddLambda = self.arcLengthController.computeDDLambda(dU, ddU_0, ddU_f, timeStep, self.theDofManager)

            # assemble total solution
            ddU = ddU_0 + ddLambda * ddU_f

            dU += ddU
            dLambda += ddLambda

            iterationCounter += 1

        self.Lambda += dLambda
        self.dLambda = dLambda

        self.arcLengthController.finishIncrement(U_n, dU, dLambda, timeStep, self.theDofManager)
        model.additionalParameters["additionalParameters"] = self.Lambda

        for stepActionType in stepActions.values():
            for action in stepActionType.values():
                action.applyAtIncrementStart(model, timeStep)

        self.journal.message(
            "Current load parameter: lambda={:6.2e}, last increment: dLambda={:6.2e}".format(self.Lambda, self.dLambda),
            self.identification,
            level=1,
        )

        return U_np, dU, P, iterationCounter, incrementResidualHistory

    def extrapolateLastIncrement(
        self,
        extrapolation: str,
        timeStep: TimeStep,
        dU: DofVector,
        dirichlets: list,
        prevTimeStep: TimeStep,
        model: FEModel,
        dLambda: float = None,
    ) -> tuple[DofVector, bool]:
        """Depending on the current setting, extrapolate the solution of the last increment.
        Also accounts for the arc length loading parameter.

        Parameters
        ----------
        extrapolation
            The type of extrapolation.
        timeStep
            The current time increment.
        dU
            The last solution increment.
        dirichlets
            The list of active dirichlet boundary conditions.
        prevTimeStep
            The previous TimeStep

        Returns
        -------
        tuple[DofVector,bool]
            - The extrapolated solution increment.
            - True if an extrapolation was performed.
        """

        # incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment

        if dLambda is None:
            # arclength control is inactive
            return super().extrapolateLastIncrement(extrapolation, timeStep, dU, dirichlets, prevTimeStep, model)

        elif extrapolation == "linear" and prevTimeStep and prevTimeStep.stepProgressIncrement:
            dLambda = dLambda * (timeStep.stepProgressIncrement / prevTimeStep.stepProgressIncrement)
            dU, isExtrapolatedIncrement = super().extrapolateLastIncrement(
                extrapolation, timeStep, dU, {}, prevTimeStep, model
            )
        else:
            dLambda = 0.0
            dU[:] = 0.0
            isExtrapolatedIncrement = False

        return dU, isExtrapolatedIncrement, dLambda

    def applyStepActionsAtStepEnd(self, model, stepActions: dict[str, StepActionBase]):
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

        if self.arcLengthController is not None:
            for stepAction in stepActions["nodeforces"].values():
                stepAction.applyAtStepEnd(model, stepMagnitude=self.Lambda)
            for stepAction in stepActions["distributedload"].values():
                stepAction.applyAtStepEnd(model, stepMagnitude=self.Lambda)
            for stepAction in stepActions["bodeforce"].values():
                stepAction.applyAtStepEnd(model, stepMagnitude=self.Lambda)

        return super().applyStepActionsAtStepEnd(model, stepActions)
