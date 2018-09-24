#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""
import numpy as np
from collections import OrderedDict, defaultdict
from fe.config.phenomena import getFieldSize, domainMapping
from fe.config.generators import getGeneratorByName
from fe.config.stepactions import stepActionFactory
from fe.config.outputmanagers import getOutputManagerByName
from fe.config.solvers import getSolverByName
from fe.utils.fieldoutput import FieldOutputController
from fe.utils.misc import  filterByJobName, stringDict
from fe.utils.plotter import Plotter
from fe.utils.exceptions import StepFailed
from fe.config.configurator import loadConfiguration, updateConfiguration
from fe.journal.journal import Journal
from fe.utils.caseinsensitivedict import CaseInsensitiveDict
from fe.utils.abqmodelconstructor import AbqModelConstructor
from time import time as getCurrentTime

def assignFieldDofIndices(nodes, constraints, domainSize):
    """ Loop over all nodes to generate the global field-dof indices.
    output is a tuple of:
        - number of total DOFS
        - orderedDict( (mechanical, indices), 
                       (nonlocalDamage, indices)
                       (thermal, indices)
                       ...)."""
        
    fieldIdxBase = 0
    fieldIndices = OrderedDict()
    for node in nodes.values():
            #delete all fields of a node, which are not active
        for field, enabled in list(node.fields.items()):
            if not enabled:
                del node.fields[field]
                
        for field in node.fields.keys():
            fieldSize = getFieldSize(field, domainSize)
            node.fields[field] = [i + fieldIdxBase for i in range(fieldSize)]
            indexList = fieldIndices.setdefault(field, []) 
            indexList += (node.fields[field])
            fieldIdxBase += fieldSize             
                
    for field, indexList in fieldIndices.items():
        fieldIndices[field] = np.array(indexList)
        
    for constraint in constraints.values():
        # some constraints may need additional Degrees of freedom (e.g. lagrangian multipliers)
        # we create them here, and assign them directly to the constraints 
        # (In contrast to true field indices, which are not directly 
        # assigned to elements/constraints but to the nodes)
        nNeededDofs = constraint.getNumberOfAdditionalNeededDofs()
        indicesOfConstraintAdditionalDofs = [i + fieldIdxBase for i in range(nNeededDofs)  ]
        constraint.assignAdditionalGlobalDofIndices ( indicesOfConstraintAdditionalDofs )
        fieldIdxBase += nNeededDofs
        
    return fieldIdxBase, fieldIndices
    
def collectStepActionsAndOptions(step, jobInfo, modelInfo, time, U, P,  stepActions, stepOptions, journal):
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
            else:
                stepActions[module][moduleName] = stepActionFactory(module)(moduleName, options, jobInfo, modelInfo, journal)
                                                               
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
        
    # the standard 'Abaqus like' model generator is invoked unconditionally and it has direct access to the inputfile
    abqModelConstructor = AbqModelConstructor(journal)
    modelInfo = abqModelConstructor.createGeometryFromInputFile(modelInfo, inputfile)
    modelInfo = abqModelConstructor.assignSectionsFromInputFile(modelInfo, inputfile)
    
    modelInfo = abqModelConstructor.createNodeFields(modelInfo, inputfile)
    modelInfo = abqModelConstructor.createConstraintsFromInputFile(modelInfo, inputfile)
    
    # create total number of dofs and orderedDict of fieldType and associated numbered dofs
    numberOfDofs, fieldIndices = assignFieldDofIndices(modelInfo['nodes'], modelInfo['constraints'], domainSize)
    
    journal.message("total size of eq. system: {:}".format(numberOfDofs), identification, 0)
    journal.printSeperationLine()

    jobName = job.get('name', '')
    time = job.get('startTime', 0.0)
    # store job info
    jobInfo = {'domainSize':        domainSize,
               'numberOfDofs':      numberOfDofs,
               'fieldIndices':      fieldIndices,
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
    
    try:
        for step in jobSteps:
            try:
                # collect all step actions in a dictionary with key of actionType
                # and concerned values 
                stepActions, stepOptions = collectStepActionsAndOptions(step, jobInfo, modelInfo, time,  U, P, stepActions, stepOptions, journal)
                
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
        if not suppressPlots:
            plotter.show()
        return success, U, P, fieldOutputController