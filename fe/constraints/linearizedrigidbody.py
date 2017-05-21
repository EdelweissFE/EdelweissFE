#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 11:34:35 2017

@author: matthias
"""

import numpy as np

from fe.utils.misc import stringDict

class Constraint:
    """ linearized Rigid Body constraint """
    
    # linear constraints are totally independent of the solution vector,
    # and, thus, need only be evaluated once (per step)
    linearConstraint = True
    
    def __init__(self, name, definitionLines, modelInfo):
        
        self.name = name
        definition = stringDict( [ e for line in definitionLines for e in line  ] )
        
        rbNset = definition.get('nSet')
        nodeSets = modelInfo['nodeSets']
        
        
        self.nodes = nodeSets[ rbNset ]
        
        
        # 1th node becomes the rp:
            
        self.rp = self.nodes[0]
        
        self.fieldsOfNodes = [ ['displacement'] for node in self.nodes ]
        
        self.fieldsOfNodes[0] += ['rotation']
        
#        nodesNonRP = [ n for n  in self.nodes if n is not  self.rp ]
#        print(nodesNonRP )1

        nDim = modelInfo['domainSize']
        
        nConstraints = len(self.nodes) + 1 # 
        nAffectedDofs =  len(self.nodes) * nDim + 1   
        

        distances = [ n.coordinates - self.rp.coordinates for n in self.nodes ]
        distanceMagnitudes = [ d @ d for d in distances  ]

        dG_dU = np.zeros( (nConstraints, nAffectedDofs )  )
        
#        for i in range(nConstraints):
#            for j, n, d, h in enumerate(zip(self.nodes, distances, distanceMagnitudes)):
#                dG_dU[i,j] = 
            
        
#        self.constraintFunctions = []
        self.additionalDofs = nConstraints
        self.additionalDofIndices =  []
        
        self.referencePoint = None
        
        self.nDof = 2 * len(self.fieldsOfNodes) + 1 + self.additionalDofs
        self.sizeStiffness =  self.nDof * self.nDof
        
        self.stiffnessContribution = np.zeros( (self.nDof, self.nDof) )
        self.effortContribution = np.zeros( self.nDof )

    def assignAdditionalDofIndices(self, dofIndices):
        self.additionalDofIndices = dofIndices
        print("dof indices for constraint:")
        print(dofIndices)
    def getAdditionalDofIndices(self,):
        return self.additionalDofIndices
    
    
    def generateConstraintStiffness(self):
        pass
    def generateConstraintForce(self):
        pass

        
    
    