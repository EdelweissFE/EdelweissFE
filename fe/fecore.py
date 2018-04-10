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
from fe.config.stepactions import stepActionFactory
from fe.config.outputmanagers import getOutputManagerByName
from fe.config.solvers import getSolverByName
from fe.config.constraints import getConstraintByName
from fe.utils.fieldoutput import FieldOutputController
from fe.utils.misc import isInteger, filterByJobName, stringDict
from fe.utils.plotter import Plotter
from fe.utils.exceptions import StepFailed
from fe.config.configurator import loadConfiguration, updateConfiguration
from fe.journal.journal import Journal
from fe.utils.caseinsensitivedict import CaseInsensitiveDict

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
        
        
        # TO DO _------------------ BETTER PROGRAMMING
        
        
        #decide if entries are labels or existing nodeSets:
        if isInteger(elSetDefinition['data'][0][0]):
            elNumbers = [int(num) for line in elSetDefinition['data'] for num in line]

            if elSetDefinition.get('generate', False):
                generateDef = elNumbers[0:3]
                els = [elements[n] for n in np.arange(generateDef[0], generateDef[1]+1, generateDef[2], dtype=int) ]
            
            elif elSetDefinition.get('boolean', False):
                booleanDef = elSetDefinition.get('boolean')
                if booleanDef=='difference':
                    els = [n for n in elementSets[name] if n.elNumber not in elNumbers ]
                
                elif booleanDef=='union':
                    els = [n for n in elementSets[name]]
                    els += [ elements[n] for n in elNumbers]

                elif booleanDef=='intersection':
                    elNumbersBase = [n.elNumber for n in elementSets[name]]
                    els = [ elements[n] for n in list(set(elNumbers).intersection(elNumbersBase))]
                    for el in els:
                        print(el.elNumber)
                    
                else:
                    raise Exception("Undefined boolean operation!")
                
                if  elSetDefinition.get('newElSet') != name:
                    name =  elSetDefinition.get('newElSet')
                else:
                    del elementSets[name]
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
    
    # generate surfaces sets
    for surfaceDef in inputfile['*surface']:
        name = surfaceDef['name']
        sType = surfaceDef.get('type', 'element')
        surface  = {} #.get(name, {})
        if sType == 'element':
            for l in surfaceDef['data']:
                elSet, faceNumber = l
                faceNumber = int(faceNumber)
                elements = modelInfo['elementSets'][elSet]
                elements +=  surface.setdefault(faceNumber, [])
                surface[faceNumber] = elements
                
        modelInfo['surfaces'][name] = surface
        
    for constraintDef in inputfile['*constraint']:
        name = constraintDef['name']
        constraint = constraintDef['type']
        data = constraintDef['data']
        
        constraint = getConstraintByName(constraint)(name, data, modelInfo)
        
        for node, nodeFields in zip(constraint.nodes, constraint.fieldsOfNodes):
            node.fields.update( [ (f, True) for f in nodeFields]  )
        
        modelInfo['constraints'][name] = constraint
        
    return modelInfo
    
def assignSections(inputfile, elementSets):
    """ Assign properties and section properties to all elements by
    the given section definitions."""
    
    for secDef in inputfile['*section']:
        if secDef['type'] == "planeUelUmat" or secDef['type'] == "solidUelUmat":
            material = [mat for mat in inputfile['*material'] if mat['id'] == secDef['material']][0]
            if secDef['type'] == "planeUelUmat":
                uelProperties = np.asarray( [ secDef['thickness'] ], dtype=float)
            else:
                uelProperties = np.array([], dtype=float)
                
            umatProperties = np.hstack(material['data'])
            for line in secDef['data']: 
               for elSet in line: 
                   for el in elementSets[elSet]:
                       el.setProperties(uelProperties, 
                                        material['name'], 
                                        material['statevars'],
                                        umatProperties)
        else:
            raise Exception("Undefined section")

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
                
    # compact storage of the model
    for generatorDefinition in inputfile['*modelGenerator']:
        gen = generatorDefinition['generator']
        modelInfo = getGeneratorByName(gen)(generatorDefinition, modelInfo, journal)
        
    modelInfo = collectNodesAndElementsFromInput(inputfile, modelInfo)
    
    # create total number of dofs and orderedDict of fieldType and associated numbered dofs
    numberOfDofs, fieldIndices = assignFieldDofIndices(modelInfo['nodes'], modelInfo['constraints'], domainSize)
    
    modelInfo['nodeSets']['all'] = list( modelInfo['nodes'].values() )
    modelInfo['elementSets']['all'] = list ( modelInfo['elements'].values() )
    
    assignSections(inputfile, modelInfo['elementSets'])
    
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
        
