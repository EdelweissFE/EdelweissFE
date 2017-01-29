#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 13:03:09 2017

@author: matthias
"""
from fe.utils.misc import stringDict
import numpy as np

def generateDirichlet(actionDefinitionLines, jobInfo, modelInfo, time, 
                                                               stepActions, 
                                                               U, P):
    dirichletIndices = []
    dirichletDelta = []
    
    nodeSets = modelInfo['nodeSets']
    
    for dirichletLine in actionDefinitionLines:
        action = stringDict(dirichletLine)        
        field = action['field']
        for x, direction  in enumerate(['1', '2', '3']):
            if direction in action:
                directionIndices = [node.fields[field][x] for node in nodeSets[action['nSet']]]
                dirichletIndices += directionIndices
                dirichletDelta += [float(action[direction])] * len(directionIndices)
                        
    dirichlet = {}
    dirichlet['indices'] =    np.array(dirichletIndices)
    dirichlet['delta'] =      np.array(dirichletDelta)
    
    return dirichlet