#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""

import numpy as np
from collections import OrderedDict, defaultdict
from fe.elements.node import Node
from fe.config.elementlibrary import elementlibrary
from fe.config.phenomena import getFieldSize, domainMapping
from fe.config.stepactions import stepActionModules
from fe.config.outputmanagers import outputManagersLibrary
from fe.config.solvers import solverLibrary
from fe.journal.journal import Journal
from time import process_time


def collectNodesAndElementsFromInput(inputfile, domainSize):
    """ Collects nodes, elements, node sets and element sets from
    the input file. """
    
    # returns an OrderedDict of {node label: node} 
    nodeDefinitions = OrderedDict()
    for nodeDefs in inputfile['*node']:
        for defLine in nodeDefs['data']:
            label = int(defLine[0])
            coordinates = np.zeros(domainSize)
            coordinates[:] = defLine[1:]
            nodeDefinitions[label] = Node(label, coordinates, )

    # returns an OrderedDict of {element Label: element}   
    elements = OrderedDict()
    for elDefs in inputfile['*element']:
        elementType = elDefs['type']
        ElementClass = elementlibrary[elementType]

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
            
    #delete all fields of a node, which are not active
    for node in nodeDefinitions.values():
        for field, enabled in list(node.fields.items()):
            if not enabled:
                del node.fields[field]

    # generate dictionary of elementObjects belonging to a specified elementset
    # or generate elementset by generate definition in inputfile
    elementSets = {}
    elementSets['all'] = [elements[e] for e in elements]
    for elSetDefinition in inputfile['*elSet']:
        name = elSetDefinition['elSet']
        if elSetDefinition.get('generate', False):
            generateDef = elSetDefinition['data'][0][0:3]
            els = [elements[n] for n in np.arange(generateDef[0], generateDef[1]+1, generateDef[2], dtype=int) ]
        else:
            els = [elements[elNum] for elNumbers in elSetDefinition['data'] for elNum in elNumbers]
        elementSets[name] = els

    # generate dictionary of nodeObjects belonging to a specified nodeset
    # or generate nodeset by generate definition in inputfile
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
    overwritten or extended. Returns a dictionary with keys defined in 
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
    external modules."""
    
    identification ="feCore"
    
    # create job dictionary from input file like e.g.
    # job  = {'data': [], 'inputFile': 'testBeam.inp', 'domain': '2d', 'name': 'testjob'}
    job = inputfile['*job'][0]

    modelDomain = job['domain']
    domainSize = domainMapping[modelDomain]

    nodes, elements, nodeSets, elementSets = collectNodesAndElementsFromInput(inputfile, domainSize)

    # assign element properties to section
    assignSections(inputfile, elementSets)
    
    # create total number of dofs and orderedDict of fieldType and associated numbered dofs
    numberOfDofs, fieldIndices = assignFieldDofIndices(nodes, domainSize)
    
    # instance of journal class with default supressFromLevel 3
    journal = Journal(verbose = verbose)
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

    # compact storage of the model
    modelInfo = {'nodes' : nodes,
                 'elements': elements,
                 'nodeSets': nodeSets,
                 'elementSets': elementSets,}
    
    jobName = job.get('name', 'defaultJob')
    # collect all job steps in a list of stepDictionaries
    jobSteps = [step for step in inputfile['*step'] if step.get('jobName', 'defaultJob') == jobName ]
                
    # collect all output managers in a list of objects     
    outputmanagers = []
    for outputDef in [output for output in inputfile['*output']     
                                if output.get('jobName', 'defaultJob') == jobName ]:
        OutputManager = outputManagersLibrary.get(outputDef['type'].lower(), None)
        if OutputManager is not None:
            outputmanagers.append( OutputManager(outputDef['data'], jobInfo, modelInfo, journal))
    
    # generate an instance of the desired solver
    solver = solverLibrary[job.get('solver', 'NIST')](jobInfo, modelInfo, journal, outputmanagers)
    U, P = solver.initialize()
    stepActions = {}
    try:
        for step in jobSteps:
            # collect all step actions in a dictionary with key of actionType
            # and concerned values 
            stepActions = collectStepActions(step, jobInfo, modelInfo, time, stepActions, U, P)
            # solve the step 
            tic =  process_time()
            success, U, P, time = solver.solveStep(step, time, stepActions, U, P,)
            toc = process_time()
            
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
        
