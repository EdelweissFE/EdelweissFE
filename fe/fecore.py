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
Created on Tue Jan  17 19:10:42 2017

@author: Matthias Neuner
"""
from collections import OrderedDict, defaultdict
from fe.config.phenomena import domainMapping
from fe.config.generators import getGeneratorByName
from fe.config.stepactions import stepActionFactory
from fe.config.outputmanagers import getOutputManagerByName
from fe.config.solvers import getSolverByName
from fe.utils.fieldoutput import FieldOutputController
from fe.utils.misc import  filterByJobName, stringDict
from fe.utils.plotter import Plotter
from fe.utils.exceptions import StepFailed
from fe.utils.dofmanager import DofManager
from fe.config.configurator import loadConfiguration, updateConfiguration
from fe.journal.journal import Journal
from fe.utils.caseinsensitivedict import CaseInsensitiveDict
from fe.utils.abqmodelconstructor import AbqModelConstructor
from time import time as getCurrentTime
    
def collectStepActionsAndOptions(step, jobInfo, modelInfo, time, U, P,  stepActions, stepOptions, fieldOutputController, journal):
    """ Parses all the defined actions for the current step, 
    and calls the respective modules, which generate step-actions based on
    computed results, model info and job information.
    The step-actions are stored in a dictionary, which is handed to 
    solveStep() in the feCore main routine afterwards.
    The step action modules decide if old step-action definitions are 
    overwritten or extended. Returns a dictionary with keys as defined in 
    stepactions."""
    
    for actionType, *definition in step['data']:
        if actionType.startswith('options'):
            # line is a option 
            options = stringDict(definition)
            category = options['category']
            stepOptions[category].update(options)
        else:
            #line is an action
            module = actionType.lower()
            options = stringDict(definition)
            moduleName = options.get('name', module)
           
            if moduleName in stepActions[module]:
                stepActions[module][moduleName].updateStepAction(options)
                journal.message("Stepaction \"{:}\" will be updated".format(moduleName), "stepActionManager", 1)
            else:
                stepActions[module][moduleName] = stepActionFactory(module)(moduleName, options, jobInfo, modelInfo, fieldOutputController, journal)
                                                               
    return  stepActions, stepOptions

def finitElementSimulation(inputfile, verbose=False, suppressPlots=False):
    """ This is the core of the finite element analysis:
    It assembles the model, and controls the respective solver based
    on the defined simulation steps.
    For each step, the step-actions (dirichlet, nodeforces) are collected by
    external modules."""
    
    identification = "feCore"
    
    journal = Journal(verbose = verbose)
    
    # create job dictionary from input file like e.g.
    job = inputfile['*job'][0]

    modelDomain = job['domain']
    domainSize = domainMapping[modelDomain]
    
    modelInfo = {'nodes' :          OrderedDict(),
                 'elements':        OrderedDict(),
                 'nodeSets':        {},
                 'elementSets':     {},
                 'surfaces':        {},
                 'constraints':     {},
                 'domainSize' :     domainSize}
                
    # call individual optional model generators
    for generatorDefinition in inputfile['*modelGenerator']:
        gen = generatorDefinition['generator']
        modelInfo = getGeneratorByName(gen)(generatorDefinition, modelInfo, journal)
        
    # the standard 'Abaqus like' model generator is invoked unconditionally, and it has direct access to the inputfile
    abqModelConstructor = AbqModelConstructor(journal)
    modelInfo = abqModelConstructor.createGeometryFromInputFile(modelInfo, inputfile)
    modelInfo = abqModelConstructor.assignSectionsFromInputFile(modelInfo, inputfile)
    modelInfo = abqModelConstructor.createConstraintsFromInputFile(modelInfo, inputfile)
    
    # create total number of dofs and orderedDict of fieldType and associated numbered dofs
    dofManager = DofManager(modelInfo)
    
    journal.message("total size of eq. system: {:}".format(dofManager.nDof), identification, 0)
    journal.printSeperationLine()

    jobName = job.get('name', '')
    time = job.get('startTime', 0.0)
    # store job info
    jobInfo = {'dofManager':        dofManager,
               'computationTime':   0.0}
               
    # add or update additional job info such as inputfile, domain, name
    jobInfo.update(job)
    jobInfo = loadConfiguration(jobInfo)
    for updateConfig in inputfile['*updateConfiguration']:
        updateConfiguration(updateConfig, jobInfo)

    # collect all job steps in a list of stepDictionaries
    jobSteps = filterByJobName(inputfile['*step'], jobName)
                
    plotter = Plotter(journal, inputfile)
    
    fieldOutputController = FieldOutputController(modelInfo, inputfile, journal)
    
    # collect all output managers in a list of objects   
    outputmanagers = []
    for outputDef in filterByJobName(inputfile['*output'], jobName):
        OutputManager = getOutputManagerByName(outputDef['type'].lower())
        managerName = outputDef.get('name', jobName + outputDef['type'] )
        definitionLines = outputDef['data']
        outputmanagers.append(OutputManager(managerName, definitionLines, jobInfo, modelInfo, 
                                            fieldOutputController, journal, plotter))
    
    # generate an instance of the desired solver
    Solver =    getSolverByName(job.get('solver','NIST'))
    solver =    Solver(jobInfo, modelInfo, journal, fieldOutputController, outputmanagers)
    U, P =      solver.initialize()
    
    stepActions = defaultdict(CaseInsensitiveDict)
    stepOptions = defaultdict(CaseInsensitiveDict)

    fieldOutputController.initializeJob(time, U, P)
    
    try:
        for step in jobSteps:
            try:
                # collect all step actions in a dictionary with key of actionType
                # and concerned values 
                stepActions, stepOptions = collectStepActionsAndOptions(step, jobInfo, modelInfo, time,  U, P, stepActions, stepOptions, fieldOutputController, journal)
                
                fieldOutputController.initializeStep(step, stepActions, stepOptions)
                for manager in outputmanagers: 
                    manager.initializeStep(step, stepActions, stepOptions)
                    
                # solve the step 
                tic =  getCurrentTime()
                success, U, P, time = solver.solveStep(step, time, stepActions, stepOptions, U, P,)
                toc = getCurrentTime()
                
                stepTime = toc - tic
                jobInfo['computationTime'] += stepTime
                
                if not success:
                    raise StepFailed()
                
                journal.printTable( [ ("Step computation time", "{:10.4f}s".format(stepTime)), ], identification, level=0)
                
            finally:
                fieldOutputController.finalizeStep(U,P)
                for manager in outputmanagers:
                    manager.finalizeStep(U, P,)
                
    except KeyboardInterrupt:
        print('')
        journal.errorMessage("Interrupted by user", identification)
        
    except StepFailed:
        journal.errorMessage("Step not finished", identification)
        
    finally:
        journal.printTable( [ ("Job computation time", "{:10.4f}s".format(jobInfo['computationTime'])), ], 
                             identification, level=0,
                             printHeaderRow=False)
        # let all output managers finalize the job
        fieldOutputController.finalizeJob(U, P)
        for manager in outputmanagers:
            manager.finalizeJob(U, P,)

        plotter.finalize()
        if not suppressPlots:
            plotter.show()
        return success, U, P, fieldOutputController
