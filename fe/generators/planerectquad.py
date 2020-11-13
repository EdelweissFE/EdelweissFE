#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:41:51 2017

@author: matthias

A mesh generator, for rectangular geometries and structured quad meshes:
    
    <-----l----->
     nX elements
     __ __ __ __
    |__|__|__|__|  A
    |__|__|__|__|  |
    |__|__|__|__|  | h
    |__|__|__|__|  | nY elements
  | |__|__|__|__|  |
  | |__|__|__|__|  V
x0|_____
  y0
  
nSets, elSets, surface : 'name'_top, _bottom, _left, _right, ...
are automatically generated

Datalines:
"""

documentation = {'x0' : '(optional) origin at x axis',
                 'y0' : '(optional) origin at y axis',
                 'h' : '(optional) height of the body',
                 'l' : '(optional) length of the body',
                 'nX': '(optional) number of elements along x',
                 'nY': '(optional) number of elements along y',
                 'elType': 'type of element'}

from fe.elements.node import Node
from fe.config.elementlibrary import getElementByName
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
    elType = getElementByName(options['elType'])
    
    testEl = elType(options['elType'], [], 0, )
    if testEl.nNodes == 4:
        nNodesX = nX+1
        nNodesY = nY+1
        
    if testEl.nNodes == 8:
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
    
    elements = []
    for x in range(nX):
        for y in range(nY):
            if testEl.nNodes == 4:
                newEl =  elType(options['elType'],
                        [ nG[x,y], nG[x+1,y], nG[x+1,y+1], nG[x, y+1] ] , currentElementLabel  ) 
                    
            elif testEl.nNodes == 8:
                newEl =  elType(options['elType'],
                        [nG[2*x,2*y], nG[2*x+2,2*y], nG[2*x+2,2*y+2], nG[2*x, 2*y+2],
                                 nG[2*x+1,2*y], nG[2*x+2,2*y+1], nG[2*x+1,2*y+2], nG[2*x, 2*y+1],
                                 ] , currentElementLabel  )
            elements.append(newEl)
            modelInfo['elements'][currentElementLabel] = newEl
            
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
    
    #element sets
    elGrid = np.asarray(elements).reshape(nX, nY)
    modelInfo['elementSets']['{:}_bottom'.format(name)] =   [ e for e in elGrid[:,0] ]  
    modelInfo['elementSets']['{:}_top'.format(name)] =      [e for e in elGrid[:,-1] ] 
    modelInfo['elementSets']['{:}_central'.format(name)] =  [ elGrid[int(nX/2),int(nY/2)]  ]
    modelInfo['elementSets']['{:}_right'.format(name)] =    [e for e in elGrid[-1,:] ]  
    modelInfo['elementSets']['{:}_left'.format(name)] =     [e for e in elGrid[0,:] ]  
    
    nShearBand = min( nX, nY )
    if nShearBand > 3:
        shearBand = [ elGrid[ int ( nX/2 + i - nShearBand/2 ), int ( nY/2 + i - nShearBand/2 ) ]  for i in range(nShearBand) ]
        modelInfo['elementSets']['{:}_shearBand'.format(name)] =  [e for e in shearBand]   
        modelInfo['elementSets']['{:}_shearBandCenter'.format(name)] =  [e for e in shearBand [  int(nShearBand/2) - 1 : int(nShearBand/2) + 2  ] ]   
    
    modelInfo['elementSets']['{:}_sandwichHorizontal'.format(name)] = []
    for elList in elGrid[1:-1,:]:
        for e in elList:
            modelInfo['elementSets']['{:}_sandwichHorizontal'.format(name)].append(e)
    
    modelInfo['elementSets']['{:}_sandwichVertical'.format(name)] = []
    for elList in elGrid[:,1:-1]:
        for e in elList:
            modelInfo['elementSets']['{:}_sandwichVertical'.format(name)].append(e)
    
    modelInfo['elementSets']['{:}_core'.format(name)] = []
    for elList in elGrid[1:-1,1:-1]:
        for e in elList:
            modelInfo['elementSets']['{:}_core'.format(name)].append(e)

    #surfaces
    modelInfo['surfaces']['{:}_bottom'.format(name)] =  {1: [e for e in elGrid[:,0] ]  }
    modelInfo['surfaces']['{:}_top'.format(name)] =     {3: [e for e in elGrid[:,-1] ]  }
    modelInfo['surfaces']['{:}_right'.format(name)] =   {2: [e for e in elGrid[-1,:] ]  }
    modelInfo['surfaces']['{:}_left'.format(name)] =    {4: [e for e in elGrid[0,:] ]  }
    
    return modelInfo
        
