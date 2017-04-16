#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""
import numpy as np
from collections import OrderedDict, defaultdict
from fe.elements.node import Node
from fe.config.elementlibrary import getElementByName
from fe.config.phenomena import getFieldSize, domainMapping
from fe.config.generators import getGeneratorByName
from fe.config.stepactions import getStepActionGeneratorByName
from fe.config.outputmanagers import getOutputManagerByName
from fe.config.solvers import solverLibrary
from fe.utils.misc import isInteger, filterByJobName
from fe.config.configurator import loadConfiguration, updateConfiguration
from fe.journal.journal import Journal
from time import time as getCurrentTime

def collectNodesAndElementsFromInput(inputfile, modelInfo):
    """ Collects nodes, elements, node sets and element sets from
    the input file. """
    
    domainSize = modelInfo['domainSize']
    
    # returns an OrderedDict of {node label: node} 
    nodeDefinitions = modelInfo['nodes']
    for nodeDefs in inputfile['*node']:
        for defLine in nodeDefs['data']:
            label = int(defLine[0])
            coordinates = np.zeros(domainSize)
            coordinates[:] = defLine[1:]
            nodeDefinitions[label] = Node(label, coordinates, )

    # returns an OrderedDict of {element Label: element}   
    elements = modelInfo['elements']
    for elDefs in inputfile['*element']:
        elementType = elDefs['type']
        ElementClass = getElementByName(elementType)

        for defLine in elDefs['data']:
            label = defLine[0]
            # store nodeObjects in elNodes list
            elNodes =  [ nodeDefinitions[n] for n in defLine[1:] ]
            newEl = ElementClass(elNodes, label)
            for iNode, node in enumerate(elNodes):
                # update node.fields dictionary with available fields from phenomena, e.g
                # OrderedDict : {'mechanical': True, 'thermal': False , ... }
                node.fields.update( [ (f, True) for f in newEl.fields[iNode] ]  )
            elements[label] = newEl
            
    # generate dictionary of elementObjects belonging to a specified elementset
    # or generate elementset by generate definition in inputfile
    elementSets = modelInfo['elementSets']
    for elSetDefinition in inputfile['*elSet']:
        name = elSetDefinition['elSet']
        
        #decide if entries are labels or existing nodeSets:
        if isInteger(elSetDefinition['data'][0][0]):
            elNumbers = [int(num) for line in elSetDefinition['data'] for num in line]
            if elSetDefinition.get('generate', False):
                generateDef = elNumbers[0:3]
                els = [elements[n] for n in np.arange(generateDef[0], generateDef[1]+1, generateDef[2], dtype=int) ]
            else:
                els = [elements[elNum] for elNum  in elNumbers]
            elementSets[name] = els
        else:
            elementSets[name]  = []
            for line in elSetDefinition['data']:
                for elSet in line:
                    elementSets[name] += elementSets[elSet]

    # generate dictionary of nodeObjects belonging to a specified nodeset
    # or generate nodeset by generate definition in inputfile
    nodeSets = modelInfo['nodeSets']
    for nSetDefinition in inputfile['*nSet']:
        name = nSetDefinition['nSet']
        
        if isInteger(nSetDefinition['data'][0][0]):
            nodes = [int(n) for line in  nSetDefinition['data'] for n in line]
            if nSetDefinition.get('generate', False):
                generateDef = nodes #nSetDefinition['data'][0][0:3]
                nodes = [nodeDefinitions[n] for n in np.arange(generateDef[0], generateDef[1]+1, generateDef[2], dtype=int) ]
            else:
                nodes = [nodeDefinitions[n] for n in nodes]
            nodeSets[name] = nodes 
        else:
            nodeSets[name]  = []
            for line in nSetDefinition['data']:
                for nSet in line:
                    nodeSets[name] += nodeSets[nSet]
    
    return modelInfo
    
def assignSections(inputfile, elementSets):
    """ Assign properties and section properties to all elements by
    the given section definitions."""
    
    for secDef in inputfile['*section']:
        if secDef['type'] == "planeUelUmat":
            material = [mat for mat in inputfile['*material'] if mat['id'] == secDef['material']][0]
            uelProperties = np.append(np.asarray(secDef['thickness']) , values = np.hstack(material['data']) )
            for line in secDef['data']: 
               for elSet in line: 
                   for el in elementSets[elSet]:
                       el.setProperties(uelProperties, material['name'], material['statevars'])

def assignFieldDofIndices(nodes, domainSize):
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
        
    return fieldIdxBase, fieldIndices
    
def collectStepActions(step, jobInfo, modelInfo, time, stepActions, U, P):
    """ Parses all the defined actions for the current step, 
    and calls the respective modules, which generate step-actions based on
    computed results, model info and job information.
    The step-actions are stored in a dictionary, which is handed to 
    solveStep() in the feCore main routine afterwards.
    The step action modules decide if old step-action definitions are 
    overwritten or extended. Returns a dictionary with keys as defined in 
    stepactions."""
    
    # create a default dictionary of type list with key defined by module action
    # and values definitions in the form of list e.g.:
    # defaultdict (list) { 'dirichlet' [ dirichlet defintion 1, dirichlet definition 2, ... ]}
    actions = defaultdict(list)
    for actionDefLine in step['data']:
        moduleName = actionDefLine[0]
        actionDefinition = actionDefLine[1:]
        actions[moduleName].append(actionDefinition)
    
    for moduleName, actionDefinitionLines in actions.items():
        #generating step-actions by external modules
        generateStepAction = getStepActionGeneratorByName(moduleName)
        stepActions[moduleName] =  generateStepAction(actionDefinitionLines, 
                                                               jobInfo, 
                                                               modelInfo, 
                                                               time, 
                                                               stepActions, 
                                                               U, P)
    return  stepActions
        
    
def finitElementSimulation(inputfile, verbose=False):
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
                 'domainSize' :     domainSize}
                
    # compact storage of the model
    for generatorDefinition in inputfile['*modelGenerator']:
        gen = generatorDefinition['generator']
        modelInfo = getGeneratorByName(gen)(generatorDefinition, modelInfo, journal)
        
    modelInfo = collectNodesAndElementsFromInput(inputfile, modelInfo)
        
    # create total number of dofs and orderedDict of fieldType and associated numbered dofs
    numberOfDofs, fieldIndices = assignFieldDofIndices(modelInfo['nodes'], domainSize)
    
    modelInfo['nodeSets']['all'] = list( modelInfo['nodes'].values() )
    modelInfo['elementSets']['all'] = list ( modelInfo['elements'].values() )
    
    assignSections(inputfile, modelInfo['elementSets'])
    
    journal.message("total size of eq. system: {:}".format(numberOfDofs), identification, 0)
    journal.printSeperationLine()

    time = job.get('startTime', 0.0)
    
    # store some job info
    jobInfo = {'domainSize':        domainSize,
               'numberOfDofs':      numberOfDofs,
               'fieldIndices':      fieldIndices,
               'computationTime':   0.0}
               
    # add or update additional job info such as inputfile, domain, name
    jobInfo.update(job)
    
    jobInfo = loadConfiguration(jobInfo)
    for updateConfig in inputfile['*updateConfiguration']:
        updateConfiguration(updateConfig, jobInfo)


    jobName = job.get('name', 'defaultJob')
    # collect all job steps in a list of stepDictionaries
    jobSteps = filterByJobName(inputfile['*step'], jobName)
                
    # collect all output managers in a list of objects     
    outputmanagers = []
    for outputDef in filterByJobName(inputfile['*output'], jobName):
        OutputManager = getOutputManagerByName(outputDef['type'].lower())
        managerName = outputDef.get('name', 'defaultName')
        definitionLines = outputDef['data']
        outputmanagers.append(OutputManager(managerName, definitionLines, jobInfo, modelInfo, journal))
    
    # generate an instance of the desired solver
    Solver =    solverLibrary[job.get('solver','NIST')]
    solver =    Solver(jobInfo, modelInfo, journal, outputmanagers)
    U, P =      solver.initialize()
    stepActions = {}
    
    try:
        for step in jobSteps:
            # collect all step actions in a dictionary with key of actionType
            # and concerned values 
            stepActions = collectStepActions(step, jobInfo, modelInfo, time, stepActions, U, P)
            
            for manager in outputmanagers: 
                manager.initializeStep(step, stepActions)
                
            # solve the step 
            tic =  getCurrentTime()
            success, U, P, time = solver.solveStep(step, time, stepActions, U, P,)
            toc = getCurrentTime()
            
            stepTime = toc - tic
            jobInfo['computationTime'] += stepTime
            
            journal.message("Step finished in {:} s".format(stepTime), identification, level=0)
            
            # let all outputmanagers finalize the step
            for manager in outputmanagers:
                manager.finalizeStep()
                
            if not success:
                journal.errorMessage("Aborting job execution", identification)
                break
            
    except KeyboardInterrupt:
        journal.errorMessage("Interrupted by user", identification)
        
    finally:
        # let all output managers finalize the job
        for manager in outputmanagers:
            manager.finalizeJob()
        journal.message("Job computation time: {:} s".format(jobInfo['computationTime']), identification, level=0)
        
