#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 11:34:35 2017

@author: matthias

Linearized rigid body constraint:
    
    Constrains a nodeset to a reference point.
    Currently only available for spatialdomain = 2D.
    Keywords are defined via datalines.
"""

import numpy as np

from fe.utils.misc import stringDict
from fe.utils.exceptions import WrongDomain

documentation = {
        'nSet': '(slave) node set, which is constrained to the reference point',
        'referencePoint': '(master) reference point'
        }

class Constraint:
    """ linearized Rigid Body constraint """
    
    # linear constraints are independent of the solution vector,
    # and, thus, need only be evaluated once (per step)
    linearConstraint = True
    
    def __init__(self, name, definitionLines, modelInfo):
        
        if modelInfo['domainSize'] != 2 :
            raise WrongDomain('Liniearized Rigid Body is currently only available for 2d domain size')
            
        self.name = name
        definition = stringDict( [ e for line in definitionLines for e in line  ] )
        
        rbNset = definition['nSet']
        nodeSets = modelInfo['nodeSets']
        
        self.rp = nodeSets [ definition['referencePoint'] ] [0]
        self.slaveNodes = nodeSets[ rbNset ]  # may also contain the RP, doesn't really matter as we remove it
        
        if self.rp in self.slaveNodes: # remove the rp from the slave node set
            self.slaveNodes = [s for s in self.slaveNodes if s is not self.rp]
        
        # all nodes
        self.nodes = self.slaveNodes +  [self.rp]
        
        nSlaves = len(self.slaveNodes)
        self.slaveNodesFields = [ ['displacement'] ] * nSlaves
        self.referencePointFields = [ ['displacement', 'rotation']  ]
        self.fieldsOfNodes = self.slaveNodesFields + self.referencePointFields

        nDim = modelInfo['domainSize']
        
        nConstraints = nSlaves * 2 # one for distance, one for angle
        nAffectedDofs =  nDim * (nSlaves + 1 ) + 1   

        distances = [ s.coordinates - self.rp.coordinates for s in self.slaveNodes ]
        dMagnitudeSquares = [ d @ d for d in distances  ]
        
        self.nDof = nAffectedDofs + nConstraints
        self.sizeStiffness =  self.nDof * self.nDof
        
        dG_dU = np.zeros( (nConstraints, nAffectedDofs )  )
        
        """
               |n1x         n1y         n2x         n2y    ...  rpx         rpy         rpA|
        -------+------------------------------------------ ... ----------------------------+
        dg1_dU |dx1         dy1         0           0      ...  -dx1        -dy1        0  |
        dg2_dU |0           0           dx2         dy2    ...  -dx2        -dy2        0  |
        ...    :...........................................................................:
        dg1A_dU|-d1y/h²     d1x/h²      0           0      ...   d1y/h²     -d1x/h²     -1.|
        dg2A_dU|0           0           -d1y/h²     d1x/h² ...   d1y/h²      d1x/h²     -1.|
        ...    |...........................................................................|
        """
        
        # derivatives distance
        for i, d in enumerate (distances):
            dG_dU[i, i*nDim  :  i*nDim + nDim ] =       d.T
            dG_dU[i, - (nDim+1)  :  -1] =            -  d.T
        # derivatives angle
        for i, (d, h2) in enumerate(zip( distances, dMagnitudeSquares )):
            x = d.T/h2
            x[1] *= -1
            dG_dU[nSlaves + i, i*nDim  :  i*nDim + nDim ] =  x[::-1].T
            dG_dU[nSlaves + i, - (nDim+1)  :  -1        ] = -x[::-1].T            
        dG_dU[nSlaves:, -1] = -1
        
        
        K = np.zeros( (self.nDof, self.nDof) )
        
        """
        K =     |   0       dG_dU.T |
                |   dG_dU   0       |
        """
        
        K[0 : dG_dU.shape[1]    , -dG_dU.shape[0]: ] =    dG_dU.T
        K[   -dG_dU.shape[0]:   , 0 : dG_dU.shape[1]: ] = dG_dU
        
        self.K = K
#        self.P = np.zeros( self.nDof ) # no contribution to RHS
        
        self.additionalGlobalDofIndices =  []
        self.nConstraints = nConstraints
    
    def getNumberOfAdditionalNeededDofs(self):
        return self.nConstraints

    def assignAdditionalGlobalDofIndices(self, additionalGlobalDofIndices):
        self.additionalGlobalDofIndices = additionalGlobalDofIndices
        
        self.globalDofIndices = np.asarray([i for node, nodeFields in zip(self.nodes, self.fieldsOfNodes) 
                                for nodeField in nodeFields  
                                    for i in node.fields[nodeField]] + self.additionalGlobalDofIndices)
                    

        
    
    