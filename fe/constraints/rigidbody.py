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
        
        self.indicesOfSlaveNodesInP  = [[i*nDim + j                 for j in range ( nDim ) ] for i in range( nSlaves ) ]
        self.indicesOfRPUinP =          [ nSlaves * nDim + j        for j in range (nDim  )]
        self.indicesOfRPPhiInP =        [ nSlaves * nDim + nDim + j for j in range (nRot )]
        
        # all nodes
        self.nodes = self.slaveNodes +  [self.referencePoint]
        
        self.slaveNodesFields = [ ['displacement'] ] * nSlaves
        self.referencePointFields = [ ['displacement', 'rotation']  ]
        self.fieldsOfNodes = self.slaveNodesFields + self.referencePointFields
       
        nConstraints = nSlaves * nDim
        nRp = 1
        
        self.nDofsOnNodes = nDim * (nSlaves + nRp  )  + nRot

        self.distancesSlaveNodeRP = [ s.coordinates - self.referencePoint.coordinates for s in self.slaveNodes ]
        
        self.nDof = self.nDofsOnNodes + nConstraints
        self.sizeStiffness =  self.nDof * self.nDof
        
        self.additionalGlobalDofIndices =  []
        self.nConstraints = nConstraints
#        self.clearLambda = False
        
    def getFields(self):
        # TODO: use
        return self.slaveNodesFields + self.referencePointFields
    
#    def acceptLastState(self,):
#        self.clearLambda = True
    
    def getNumberOfAdditionalNeededDofs(self):
        return self.nConstraints

    def assignAdditionalGlobalDofIndices(self, additionalGlobalDofIndices):
        self.additionalGlobalDofIndices = additionalGlobalDofIndices
        
        self.globalDofIndices = np.asarray([i for node, nodeFields in zip(self.nodes, self.fieldsOfNodes) 
                                for nodeField in nodeFields  
                                    for i in node.fields[nodeField]] + self.additionalGlobalDofIndices)
    
    def Rphi_3D_x_(self, phi, d):
        phi = phi + np.pi / 2  * d
        i  =0.0 if d > 0 else 1.0
        return np.array([[i,            0,            0,           ],
                         [0,            np.cos(phi),  -np.sin(phi) ],
                         [0,            np.sin(phi),  +np.cos(phi) ],])

    def Rphi_3D_y_(self, phi, d):
        phi = phi + np.pi / 2  * d
        i  =0.0 if d > 0 else 1.0
        return np.array([[np.cos(phi),   0,           +np.sin(phi),],
                          [0,            i,           0            ],
                          [-np.sin(phi), 0,           +np.cos(phi) ],])
    
    def Rphi_3D_z_(self, phi, d):
        phi = phi + np.pi / 2  * d
        i  =0.0 if d > 0 else 1.0
        return np.array([[ np.cos(phi),  -np.sin(phi),   0  ],
                         [ np.sin(phi),  +np.cos(phi),   0  ],
                         [0,             0,              i,],])

    def applyConstraint(self, Un1, PExt, V, increment):
        
        nConstraints    = self.nConstraints
        nDim            = self.nDim
        
        nU = self.nDof - nConstraints
        nSlaves = len ( self.slaveNodes )
        

        URp =       Un1[self.indicesOfRPUinP]
        PhiRp =     Un1[self.indicesOfRPPhiInP]
        Lambdas =   Un1[nU:].reshape( (nDim, -1), order='F')
        
#        if self.clearLambda:
#            Lambdas[:] = 0.0
#            self.clearLambda = False
        
        K = V.reshape( self.nDof, self.nDof, order='F')

        G = np.zeros( (nDim, 9) )
        H = np.zeros( (nDim, 9, 9) )
        
        # d/dUN
        G[:, 0:3 ] = - np.identity(nDim) 
        G[:, 3:6 ] = + np.identity(nDim) 
        
        Rx, Ry, Rz = self.Rphi_3D_x_, self.Rphi_3D_y_, self.Rphi_3D_z_
        
        RotationMatricesAndDerivatives = [ [R(phi, d)  for d in range(3) ] for R, phi in zip ( (Rx,Ry,Rz), PhiRp)  ]
        Rx = RotationMatricesAndDerivatives[0]
        Ry = RotationMatricesAndDerivatives[1]
        Rz = RotationMatricesAndDerivatives[2]
        
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

        
        #start and end of Lambda in P
        L0, LF = nU, nU + nDim
        
        for i in range( nSlaves ):
        
            d0 = self.distancesSlaveNodeRP[i]
#            print(d0)
            Lambda = Lambdas[:,i]
            indcsUNode = self.indicesOfSlaveNodesInP[i]
            
            Un =     Un1[indcsUNode]

            g =  -d0 - (Un - URp) + Rz[0] @ Ry[0] @ Rx[0] @ d0
            
            # TODO: maybe implemented more performant using np.tensordot (but also morge ugly..)            
            for i in range(3):
                G[:, 2*nDim +i] = RDerivativeProductsI[i] @ d0
                for j in range(3):
                    H[:, 6+i, 6+j] = RDerivativeProductsII[i][j] @ d0 
            indcsU = np.array( indcsUNode + self.indicesOfRPUinP + self.indicesOfRPPhiInP, dtype=np.int)
            
            KUU = K[0:nU, 0:nU]
            KUL = K[0:nU,   L0 : LF ]
            KLU = K[L0: LF , 0:nU ]

            PExt[indcsU]    -= Lambda.T @ G
            PExt[L0:LF]     -= g 
            
            #for coordinate like - black magic - access in KUU..:
            t = np.tile( indcsU , (indcsU.shape[0] ,1))
            
            KUL [ indcsU, : ] += G.T
            KLU [ :, indcsU ] += G
            KUU[ t.T, t ]     += np.tensordot( Lambda, H, (0,0) ) # L_[i] H_[i,j,k]
            # ... corresponds to:
#            for i in range(6,9):
#                for j in range(6,9):
#                    for k in range(nDim) :
#                        KUU [ indcsU[i], indcsU[j] ]  += Lambda[k] * H[k,i,j]
            

            L0 += nDim
            LF += nDim
#            print(Lambda)
