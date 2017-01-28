#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""

import numpy as np
from collections import OrderedDict, defaultdict
from fe.config.elementlibrary import elementlibrary

from fe.elements.node import Node
from fe.config.phenomena import getFieldSize, domainMapping
from fe.config.stepactions import stepActionModules
from fe.config.outputmanagers import outputManagersLibrary
from fe.solvers.nonlinearImplicitStatic import NIST
from fe.journal.journal import Journal
from time import process_time


def collectNodesAndElementsFromInput(inputfile, domainSize):
    """ Collects nodes, elements, nodeSets and elementSets from
    the input file. """
    
    nodeDefinitions = OrderedDict()
    for nodeDefs in inputfile['*node']:
        for defLine in nodeDefs['data']:
            label = int(defLine[0])
            coordinates = np.zeros(domainSize)
            coordinates[:] = defLine[1:]
            nodeDefinitions[label] = Node(label, coordinates, )
            
    elements = OrderedDict()
    for elDefs in inputfile['*element']:
        elementType = elDefs['type']
        ElementClass = elementlibrary[elementType]

        for defLine in elDefs['data']:
            label = defLine[0]
            elNodes =  [ nodeDefinitions[n] for n in defLine[1:] ]
            newEl = ElementClass(elNodes, label) 
            for iNode, node in enumerate(elNodes):
                node.fields.update( [ (f, True) for f in newEl.fields[iNode] ]  )
            elements[label] = newEl
            
    #disable all fields of a node, which are not active
    for node in nodeDefinitions.values():
        for field, enabled in list(node.fields.items()):
            if not enabled:
                del node.fields[field]

    elementSets = {}
    elementSets['all'] = [elements[e] for e in elements]
    for elSetDefinition in inputfile['*elSet']:
        name = elSetDefinition['elSet']
        els = [elements[e] for l in elSetDefinition['data'] for e in l]
        elementSets[name] = els
        
    nodeSets = {}
    nodeSets['all'] = [nodeDefinitions[n] for n in nodeDefinitions]
    for nSetDefinition in inputfile['*nSet']:
        name = nSetDefinition['nSet']
        if nSetDefinition.get('generate', False):
            generateDef = nSetDefinition['data'][0][0:3]
            nodes = [nodeDefinitions[n] for n in np.arange(generateDef[0], generateDef[1]+1, generateDef[2], dtype=int) ]
        else:
            nodes = [nodeDefinitions[n] for l in nSetDefinition['data'] for n in l]
        nodeSets[name] = nodes 
    
    return nodeDefinitions, elements, nodeSets, elementSets
    
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
    """ Loop over all nodes,
    to generate the global field-dof indices."""
    
    fieldIdxBase = 0
    fieldIndices = OrderedDict()
    for node in nodes.values():
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
        The step-actions are stored in a dictionary, which is handled to 
        solveStep() in the feCore main routine afterwords.
        The step action modules decide if old step-action definitions are 
        overwritten or extended."""
    
    actions = defaultdict(list)
    for actionDefLine in step['data']:
        #parsing and storing 
        moduleName = actionDefLine[0]
        actionDefinition = actionDefLine[1:]
        actions[moduleName].append(actionDefinition)
    
    for moduleName, actionDefinitionLines in actions.items():
        #generating step-actions by external modules
        actionModule = stepActionModules[moduleName]
        stepActions[moduleName] =  actionModule(actionDefinitionLines, 
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
        corresponing external modules."""
    
    identificiation ="feCore"
    
    job = inputfile['*job'][0]
    modelDomain = job['domain']
    domainSize = domainMapping[modelDomain]

    nodes, elements, nodeSets, elementSets = collectNodesAndElementsFromInput(inputfile, domainSize)    
    assignSections(inputfile, elementSets)
    numberOfDofs, fieldIndices = assignFieldDofIndices(nodes, domainSize)
    
    journal = Journal()
    
    journal.message("total size of eq. system: {:}".format(numberOfDofs), identificiation, 0)
    journal.printSeperationLine()

    time = job.get('startTime', 0.0)
    
    # store some job info
    jobInfo = {'domainSize':domainSize,
               'numberOfDofs':numberOfDofs,
               'fieldIndices':fieldIndices,
               'computationTime':0.0}
    jobInfo.update(job)
    
    # compact storage of the model
    modelInfo = {'nodes' : nodes,
                 'elements': elements,
                 'nodeSets': nodeSets,
                 'elementSets': elementSets,}
    
    # collect all job steps
    jobName = job.get('name', 'defaultJob')
    jobSteps = [step for step in inputfile['*step'] if step.get('jobName', 'defaultJob') == jobName ]
    
    # collect all output managers      
    outputmanagers = []
    for outputDef in [output for output in inputfile['*output']     
                                if output.get('jobName', 'defaultJob') == jobName ]:
        OutputManager = outputManagersLibrary.get(outputDef['type'].lower(), None)
        if OutputManager is not None:
            outputmanagers.append( OutputManager(outputDef['data'], jobInfo, modelInfo, journal))
        
    
    # generate an instance for the desired solver
    solver = NIST(jobInfo, modelInfo, journal, outputmanagers)
    U, P = solver.initialize()
    stepActions = {}
    try:
        for step in jobSteps:
            
            # collect all step actions
            stepActions = collectStepActions(step, jobInfo, modelInfo, time, stepActions, U, P)
            
            # solve the step 
            tic =  process_time()
            success, U, P, time = solver.solveStep(step, time, stepActions, U, P,)
            toc = process_time()
            
            stepTime = toc - tic
            jobInfo['computationTime'] += stepTime
            
            journal.message("Step finished in {:} s".format(stepTime), identificiation, level=0)
            
            # let all outputmanagers finalize the step
            for manager in outputmanagers:
                manager.finalizeStep()
                
            if not success:
                journal.errorMessage("Aborting job execution", identificiation)
                break
            
    except KeyboardInterrupt:
        journal.errorMessage("Interrupted by user", identificiation)
        
    finally:
        # let all output managers finalize the job
        for manager in outputmanagers:
            manager.finalizeJob()
        journal.message("Job computation time: {:} s".format(jobInfo['computationTime']), identificiation, level=0)
        
