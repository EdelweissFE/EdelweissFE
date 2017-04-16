#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 19:33:06 2017

@author: matthias
"""

from fe.utils.misc import stringDict
import numpy as np

def generateAction(actionDefinitionLines, jobInfo, modelInfo, 
                       time, stepActions, U, P):
    """ create nodeForces dictionary with nodeForce in 
        keytype 'indices': array of global dof indices
                'delta':   prescribed deltaValue """
    
    nodeForceIndices = []
    nodeForceDelta = []
    
    nodeSets = modelInfo['nodeSets']
    
    for dirichletLine in actionDefinitionLines:
        action = stringDict(dirichletLine)        
        field = action['field']
        for x, direction  in enumerate(['1', '2', '3']):
            if direction in action:
                directionIndices = [node.fields[field][x] for node in nodeSets[action['nSet']]]
                nodeForceIndices += directionIndices
                nodeForceDelta += [float(action[direction])] * len(directionIndices)
                        
    nodeForces = {}
    nodeForces['indices'] =    np.array(nodeForceIndices)
    nodeForces['delta'] =      np.array(nodeForceDelta)
    
    return nodeForces