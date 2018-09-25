#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 09:44:10 2018

@author: m9
"""

from fe.elements.node import Node
from fe.config.elementlibrary import getElementByName
from fe.utils.misc import isInteger
from fe.config.constraints import getConstraintByName

import numpy as np

class AbqModelConstructor:
    
    def __init__(self, journal):
        self.journal = journal

    def createGeometryFromInputFile(self, modelInfo, inputFile):
        """ Collects nodes, elements, node sets and element sets from
        the input file. """
        
        domainSize = modelInfo['domainSize']
        
        # returns an OrderedDict of {node label: node} 
        nodeDefinitions = modelInfo['nodes']
        for nodeDefs in inputFile['*node']:
            for defLine in nodeDefs['data']:
                label = int(defLine[0])
                coordinates = np.zeros(domainSize)
                coordinates[:] = defLine[1:]
                nodeDefinitions[label] = Node(label, coordinates, )
    
        # returns an OrderedDict of {element Label: element}   
        elements = modelInfo['elements']
        
        for elDefs in inputFile['*element']:
            elementType = elDefs['type']
            ElementClass = getElementByName(elementType)
    
            for defLine in elDefs['data']:
                label = defLine[0]
                # store nodeObjects in elNodes list
                elNodes =  [ nodeDefinitions[n] for n in defLine[1:] ]
                newEl = ElementClass(elementType, elNodes, label)
#                for iNode, node in enumerate(elNodes):
#                    # update node.fields dictionary with available fields from phenomena, e.g
#                    # OrderedDict : {'mechanical': True, 'thermal': False , ... }
#                    node.fields.update( [ (f, True) for f in newEl.fields[iNode] ]  )
                elements[label] = newEl
                
        # generate dictionary of elementObjects belonging to a specified elementset
        # or generate elementset by generate definition in inputfile
        elementSets = modelInfo['elementSets']
        
        for elSetDefinition in inputFile['*elSet']:
            name = elSetDefinition['elSet']
            
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
        for nSetDefinition in inputFile['*nSet']:
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
        
        modelInfo['nodeSets']['all'] = list( modelInfo['nodes'].values() )
        modelInfo['elementSets']['all'] = list ( modelInfo['elements'].values() )
        
        # generate surfaces sets
        for surfaceDef in inputFile['*surface']:
            name = surfaceDef['name']
            sType = surfaceDef.get('type', 'element').lower()
            surface  = {} 
            if sType == 'element':
                for l in surfaceDef['data']:
                    elSet, faceNumber = l
                    faceNumber = int(faceNumber)
                    elements = modelInfo['elementSets'][elSet]
                    elements +=  surface.setdefault(faceNumber, [])
                    surface[faceNumber] = elements
                    
            modelInfo['surfaces'][name] = surface
            
        return modelInfo

    def createConstraintsFromInputFile(self, modelInfo, inputFile): 
        for constraintDef in inputFile['*constraint']:
            name = constraintDef['name']
            constraint = constraintDef['type']
            data = constraintDef['data']
            
            constraint = getConstraintByName(constraint)(name, data, modelInfo)
            
            for node, nodeFields in zip(constraint.nodes, constraint.fieldsOfNodes):
                node.fields.update( [ (f, True) for f in nodeFields]  )
            
            modelInfo['constraints'][name] = constraint
            
        return modelInfo
    
    def assignSectionsFromInputFile(self, modelInfo, inputFile):
        """ Assign properties and section properties to all elements by
        the given section definitions."""
        
        elementSets = modelInfo['elementSets']
    
        for secDef in inputFile['*section']:
            if secDef['type'] == "planeUelUmat" or secDef['type'] == "solidUelUmat":
                material = [mat for mat in inputFile['*material'] if mat['id'] == secDef['material']][0]
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
                                            umatProperties)
            else:
                raise Exception("Undefined section")
                
        return modelInfo
    
    def createNodeFields(self, modelInfo, inputFile):
        for element in modelInfo['elements'].values():
            for node, nodeFields in zip ( element.nodes, element.fields ):
#            for iNode, node in enumerate(elNodes):
                # update node.fields dictionary with available fields from phenomena, e.g
                # OrderedDict : {'mechanical': True, 'thermal': False , ... }
#                node.fields.update( [ (f, True) for f in newEl.fields[iNode] ]  )
                node.fields.update( [ (f, True) for f in nodeFields ]  )
        return modelInfo
        