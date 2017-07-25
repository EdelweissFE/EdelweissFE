#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 21:03:23 2017

@author: matthias
"""

from collections import OrderedDict, defaultdict

def extractNodesFromElementSet(elementSet):
    nodeCounter = 0
    partNodes = OrderedDict() # node -> index in nodelist
    for element in elementSet:
        for node in element.nodes:
            # if the node is already in the dict, get its index, 
            # else insert it, and get the current idx = counter. increase the counter
            idx = partNodes.setdefault(node, nodeCounter)
            if idx == nodeCounter:
                # the node was just inserted, so increase the counter of inserted nodes
                nodeCounter+=1
    return partNodes

def disassembleElsetToEnsightShapes(elementSet):
    """
    elset -> {shape : [index in elset, ... ], }
    """
    elements = defaultdict(list)
    for i, el in enumerate(elementSet):
        elements[el.ensightType].append(i)
    return elements

def transferElsetResultsToElset( elsetTarget, elsetOrigin, resultsTarget, resultsOrigin ):
    
    for i, el in enumerate(elsetTarget):
        el.__index__in__elsetTarget = i
    indices = [ el.__index__in__elsetTarget for el in elsetOrigin]
    
    if resultsOrigin.ndim == 1:
        resultsTarget[ indices ] = resultsOrigin
    elif resultsOrigin.ndim == 2:
        resultsTarget[ indices, : ] = resultsOrigin
     
