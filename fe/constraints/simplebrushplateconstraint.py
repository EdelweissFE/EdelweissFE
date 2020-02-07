#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 09:03:42 2020

@author: matthias
"""

documentation = {
        'nSet': '(slave) node set, which is constrained to the reference point',
        'referencePoint': '(master) reference point'
        }

import numpy as np

from fe.utils.misc import stringDict

#from numba import jit

class Constraint:
    """ Simple brush plate constraint (Concrete Tests by Kupfer),
        ONLY in global principal directions!"""
    
    def __init__(self, name, definitionLines, modelInfo):
            
        self.name = name
        definition = stringDict( [ e for line in definitionLines for e in line  ] )
        
        self.nDim = modelInfo['domainSize']
        
        self.constraintDimension = int ( definition['direction'] ) - 1
        
        nDim = self.nDim
        
        rbNset =    definition['nSet']
        nodeSets =  modelInfo['nodeSets']
        
        self.referencePoint = nodeSets [ definition['referencePoint'] ] [0]
        self.slaveNodes = nodeSets[ rbNset ]  # may also contain the RP, doesn't really matter as we remove it
        
        if self.referencePoint in self.slaveNodes: # remove the rp from the slave node set
            self.slaveNodes = [s for s in self.slaveNodes if s is not self.referencePoint]
                        
        # nRot = 3
        nSlaves = len(self.slaveNodes)
        
        self.nConstraints = nSlaves
        
        dG_dU = np.zeros( (self.nConstraints, self.nConstraints + 1 )  )
        
        np.fill_diagonal(dG_dU[:, :-1] , 1.0)
        dG_dU[:, -1] = -1.0
        
        self.nDofsOnNodes = nDim * (nSlaves + 1  )  
        
        self.nDof = self.nDofsOnNodes + self.nConstraints
        
  
        self.dG_dU = dG_dU
        
        
        self.additionalGlobalDofIndices =  []
        
        # all nodes
        self.nodes = self.slaveNodes +  [self.referencePoint, ] # RP is the last one, so we find him always
        
        self.slaveNodesFields       = [ ['displacement'] ] * nSlaves
        self.referencePointFields   = [ ['displacement'],  ]
        self.fieldsOfNodes = self.slaveNodesFields + self.referencePointFields
       

    
    def getNodes(self):    
        return self.nodes
    
    def getNodeFields(self):
        raise Exception ("check what this does!")
        # TODO: use
        return self.slaveNodesFields + self.referencePointFields
    
    def getNumberOfAdditionalNeededDofs(self):
        return self.nConstraints

    def assignAdditionalGlobalDofIndices(self, additionalGlobalDofIndices):
        self.additionalGlobalDofIndices = additionalGlobalDofIndices
        
        self.globalDofIndices = np.asarray([i for node, nodeFields in zip(self.nodes, self.fieldsOfNodes) 
                                for nodeField in nodeFields  
                                    for i in node.fields[nodeField]] + self.additionalGlobalDofIndices)
        
    def assignGlobalDofIndices(self):
        pass
    
    def getGlobalDofIndices(self):
        pass
    
    def applyConstraint(self, Un1, PExt, V, increment):
        
        idcsNodeDofs = np.arange( self.constraintDimension , (self.nConstraints+1) * self.nDim, self.nDim)
        
        
        LambdaN1 = Un1[ -self.nConstraints : ]
        
        Ke = V.reshape( ( self.nDof, self.nDof) )
        
        Ke[ idcsNodeDofs, -self.nConstraints :   ] += self.dG_dU.T        
        Ke[ -self.nConstraints : ,  idcsNodeDofs ] += self.dG_dU 
        
        PExt[ idcsNodeDofs ] -= LambdaN1 .dot( self.dG_dU )
        

        
