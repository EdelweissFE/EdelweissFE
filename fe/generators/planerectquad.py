#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:41:51 2017

@author: matthias
"""

from fe.elements.node import Node
from fe.config.elementlibrary import elementlibrary
from fe.utils.misc import stringDict 

import numpy as np

def generateModelData(generatorDefinition, modelInfo, journal):
    options = generatorDefinition['data']
    options = stringDict( [ e for l in options for e in l ] )
    
    name = generatorDefinition.get('name', 'planeRect')
    
    x0 = float(options.get('x0', 0.0))
    y0 = float(options.get('y0', 0.0))
    h = float(options.get('h', 1.0))
    l = float(options.get('l', 1.0))
    nX =int(options.get('nX', 10))
    nY =int(options.get('nY', 10))
    elType = elementlibrary[options['elType']]
    
    
    if elType.nNodes == 4:
        nNodesX = nX+1
        nNodesY = nY+1
        
    if elType.nNodes == 8:
        nNodesX = 2*nX+1
        nNodesY = 2*nY+1
        
    grid =  np.mgrid[ x0: x0+l: nNodesX*1j , y0: y0+h: nNodesY*1j ,  ]
    
    nodes = []
    currentNodeLabel = 1
    
    for x in range(nNodesX):
        for y in range(nNodesY):
            node = Node(currentNodeLabel, grid[:,x,y] )
            modelInfo['nodes'][currentNodeLabel] = node
            nodes.append( node )
            currentNodeLabel +=1
    
    nG = np.asarray(nodes).reshape(nNodesX, nNodesY)
    
    currentElementLabel = 1
    
    for x in range(nX):
        for y in range(nY):
            if elType.nNodes == 4:
                newEl =  elType([ nG[x,y], nG[x+1,y], nG[x+1,y+1], nG[x, y+1] ] , currentElementLabel  ) 
            elif elType.nNodes == 8:
                newEl =  elType([nG[2*x,2*y], nG[2*x+2,2*y], nG[2*x+2,2*y+2], nG[2*x, 2*y+2],
                                 nG[2*x+1,2*y], nG[2*x+2,2*y+1], nG[2*x+1,2*y+2], nG[2*x, 2*y+1],
                                 ] , currentElementLabel  ) 
                
            modelInfo['elements'][ currentElementLabel] =newEl
            for i, node in enumerate(newEl.nodes):
                node.fields.update( [ (f, True) for f in newEl.fields[i] ]  )
            currentElementLabel +=1
        
        #nodesets:
        modelInfo['nodeSets'][ '{:}_left'.format(name) ] =  [n for n in nG[0,:]]
        modelInfo['nodeSets'][ '{:}_right'.format(name) ] = [n for n in nG[-1,:]]
        modelInfo['nodeSets'][ '{:}_top'.format(name) ] =   [n for n in nG[:,-1]]
        modelInfo['nodeSets'][ '{:}_bottom'.format(name) ] =   [n for n in nG[:,0]]
        
        modelInfo['nodeSets'][ '{:}_leftBottom'.format(name) ] =   [ nG[0,0] ]
        modelInfo['nodeSets'][ '{:}_leftTop'.format(name) ] =      [ nG[0,-1] ]
        modelInfo['nodeSets'][ '{:}_rightBottom'.format(name) ] =  [ nG[-1,0] ]
        modelInfo['nodeSets'][ '{:}_rightTop'.format(name) ] =     [ nG[-1,-1] ]
        
    return modelInfo
        
