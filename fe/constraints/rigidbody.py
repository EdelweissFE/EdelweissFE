#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 11:21:39 2018

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
    """ Rigid Body constraint """
    
    def __init__(self, name, definitionLines, modelInfo):
            
        self.name = name
        definition = stringDict( [ e for line in definitionLines for e in line  ] )
        
        self.nDim = modelInfo['domainSize']
        nDim = self.nDim
        
        rbNset =    definition['nSet']
        nodeSets =  modelInfo['nodeSets']
        
        self.referencePoint = nodeSets [ definition['referencePoint'] ] [0]
        self.slaveNodes = nodeSets[ rbNset ]  # may also contain the RP, doesn't really matter as we remove it
        
        if self.referencePoint in self.slaveNodes: # remove the rp from the slave node set
            self.slaveNodes = [s for s in self.slaveNodes if s is not self.referencePoint]
                        
            
        nRot = 3
        nSlaves = len(self.slaveNodes)
        
        self.indicesOfSlaveNodesInP  = [[       i * nDim + j        for j in range ( nDim ) ] for i in range( nSlaves ) ]
        self.indicesOfRPUinP =          [ nSlaves * nDim + j        for j in range ( nDim ) ]
        self.indicesOfRPPhiInP =        [ nSlaves * nDim + nDim + j for j in range ( nRot ) ]
        
        # all nodes
        self.nodes = self.slaveNodes +  [self.referencePoint]
        
        self.slaveNodesFields       = [ ['displacement'] ] * nSlaves
        self.referencePointFields   = [ ['displacement', 'rotation']  ]
        self.fieldsOfNodes = self.slaveNodesFields + self.referencePointFields
       
        nConstraints = nSlaves * nDim
        nRp = 1
        
        self.nDofsOnNodes = nDim * (nSlaves + nRp  )  + nRot

        self.distancesSlaveNodeRP = [ s.coordinates - self.referencePoint.coordinates for s in self.slaveNodes ]
        
        self.nDof = self.nDofsOnNodes + nConstraints
        self.sizeStiffness =  self.nDof * self.nDof
        
        self.additionalGlobalDofIndices =  []
        self.nConstraints = nConstraints
        
        self.nRot = 3
    
    def getNodes(self):    
        return self.nodes
    
    def getNodeFields(self):
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
    
    def Rz_2D(self, phi, derivative):
        phi = phi + np.pi / 2  * derivative
        return np.array([[ np.cos(phi),  -np.sin(phi) ],
                         [ np.sin(phi),  +np.cos(phi) ],])
    
    def Rx_3D(self, phi, derivative):
        phi = phi + np.pi / 2  * derivative
        i  = 0.0 if derivative > 0 else 1.0
        return np.array([[i,            0,            0,           ],
                         [0,            np.cos(phi),  -np.sin(phi) ],
                         [0,            np.sin(phi),  +np.cos(phi) ],])

    def Ry_3D(self, phi, derivative):
        phi = phi + np.pi / 2  * derivative
        i  = 0.0 if derivative > 0 else 1.0
        return np.array([[np.cos(phi),   0,           +np.sin(phi),],
                          [0,            i,           0            ],
                          [-np.sin(phi), 0,           +np.cos(phi) ],])
    
    def Rz_3D(self, phi, derivative):
        phi = phi + np.pi / 2  * derivative
        i  = 0.0 if derivative > 0 else 1.0
        return np.array([[ np.cos(phi),  -np.sin(phi),   0  ],
                         [ np.sin(phi),  +np.cos(phi),   0  ],
                         [0,             0,              i,],])

    def applyConstraint(self, Un1, PExt, V, increment):
        
        nConstraints    = self.nConstraints
        nDim            = self.nDim
        nRot            = self.nRot
        
        nU              = self.nDof - nConstraints # nDofs (disp., rot.) without Lagrangian multipliers
        nSlaves         = len ( self.slaveNodes )
        
        URp             = Un1[self.indicesOfRPUinP]
        PhiRp           = Un1[self.indicesOfRPPhiInP]
        Lambdas         = Un1[nU:].reshape( (nDim, -1), order='F')
        
        K = V.reshape( self.nDof, self.nDof, order='F')

        G = np.zeros( (nDim, 9) )
        H = np.zeros( (nDim, 9, 9) )
        
        # dg/dU_Node and # dg/dU_RP
        G[:, 0:nDim ]       = - np.identity(nDim) 
        G[:, nDim:2*nDim ]  = + np.identity(nDim) 
        
        if nDim == 3:
        
            Rx, Ry, Rz = self.Rx_3D, self.Ry_3D, self.Rz_3D
            
            RotationMatricesAndDerivatives = [ [R(phi, derivative)  for derivative in range(3) ] for R, phi in zip ( (Rx,Ry,Rz), PhiRp)  ]
            Rx = RotationMatricesAndDerivatives[0]
            Ry = RotationMatricesAndDerivatives[1]
            Rz = RotationMatricesAndDerivatives[2]
            
            T = Rz[0] @ Ry[0] @ Rx[0]
            
            RDerivativeProductsI = (
                    Rz[0] @ Ry[0] @ Rx[1],
                    Rz[0] @ Ry[1] @ Rx[0],
                    Rz[1] @ Ry[0] @ Rx[0],)
            
            RDerivativeProductsII = (
                    (Rz[0] @ Ry[0] @ Rx[2],
                     Rz[0] @ Ry[1] @ Rx[1],
                     Rz[1] @ Ry[0] @ Rx[1],),
                     
                    (Rz[0] @ Ry[1] @ Rx[1],
                     Rz[0] @ Ry[2] @ Rx[0],
                     Rz[1] @ Ry[1] @ Rx[0],),
                     
                    (Rz[1] @ Ry[0] @ Rx[1],
                     Rz[1] @ Ry[1] @ Rx[0],
                     Rz[2] @ Ry[0] @ Rx[0],),)
                    
        elif nDim == 2:
            raise Exception("rigid body constraint not yet implemented for 2D")
        
        #start and end of Lambda in P
        L0, LF = nU, nU + nDim
        
        for i in range( nSlaves ):
        
            d0          = self.distancesSlaveNodeRP[i]
            indcsUNode  = self.indicesOfSlaveNodesInP[i]
            Un          = Un1[indcsUNode]
            Lambda      = Lambdas[:,i]
            
            g =  -d0 - (Un - URp) + T @ d0
            
            for i in range(nRot):
                G[:, 2*nDim +i] = RDerivativeProductsI[i] @ d0
                for j in range(nRot):
                    # only the rot. block is nonzero
                    H[:, 2*nDim+i, 2*nDim+j] = RDerivativeProductsII[i][j] @ d0 
                    
            indcsU = np.array( indcsUNode + self.indicesOfRPUinP + self.indicesOfRPPhiInP, dtype=np.int)

            PExt[indcsU]    -= Lambda.T @ G
            PExt[L0:LF]     -= g 
            
            KUU = K[0:nU, 0:nU]
            KUL = K[0:nU, L0:LF ]
            KLU = K[L0:LF, 0:nU ]
            
            # for KUU, only the Phi_RP block is nonzero
            KUU[-nRot:, -nRot:] += np.einsum('i,ijk->jk', Lambda, H[:,-nRot:, -nRot:], ) # L_[i] H_[i,j,k]
            KUL [ indcsU, : ]   += G.T
            KLU [ :, indcsU ]   += G

            L0 += nDim
            LF += nDim
