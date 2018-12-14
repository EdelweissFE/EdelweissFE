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
from fe.utils.exceptions import WrongDomain


class Constraint:
    """ Rigid Body constraint """
    
    # linear constraints are independent of the solution vector,
    # and, thus, need only be evaluated once (per step)
    linearConstraint = True
    
    def __init__(self, name, definitionLines, modelInfo):
        
#        if modelInfo['domainSize'] != 2 :
#            raise WrongDomain('Rigid Body is currently only available for 2d domain size')
            
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

        self.nDim = modelInfo['domainSize']
        nDim = self.nDim
       
        nConstraints = nSlaves * nDim
        
        nRp = 1
        
        self.nAffectedDofs = nDim * (nSlaves + nRp + nRp ) 

        self.d0s = [ s.coordinates - self.rp.coordinates for s in self.slaveNodes ]
        
        self.nDof = self.nAffectedDofs + nConstraints
        self.sizeStiffness =  self.nDof * self.nDof
        
        self.additionalGlobalDofIndices =  []
        self.nConstraints = nConstraints
        
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
    
    
    def G_and_H_Rot(self, Unx, Uny, URpx, URpy, PhiRp, d0):
        
        """" g(...) = 0 = -PhiRp + atan2( a ,b )
        with a = d0_x * dN_y - d0y * dN_x
             b = d0_x * dN_x + d0y * dN_y
        and
            d0 = X_n - X_rp
            dN = x_n - x_rp = d0 + U_n - Urp
        """
        
        G = np.zeros ( (5) )
        H = np.zeros ( (5,5) )
        
        dN = np.copy(d0)
        dN[0] += ( Unx - URpx )
        dN[1] += ( Uny - URpy )
        
        a = d0[0] * dN[1] - d0[1] * dN[0]
        b = d0[0] * dN[0] + d0[1] * dN[1]
        
        """
        derivatives of atan2
        """
        
        b2plusa2_inv = 1. / (b**2 + a**2)
        d_atan2_d_a = b * b2plusa2_inv
        d_atan2_d_b = -a * b2plusa2_inv
        
        """
        derivatives of a, b wrt Un and Urp
        """
        
        da_dUni = [-d0[1], +d0[0], +d0[1], -d0[0]]
        db_dUni = [d0[0], d0[1], -d0[0], -d0[1]]
        
        for i in range(4):
            G[i] = d_atan2_d_a * da_dUni[i] + d_atan2_d_b * db_dUni[i]
        
        G[4] = -1
        
        """
        2nd order derivatives for H
        """
        
        b2plusa2_inv2 = b2plusa2_inv ** 2
        
        d2_atan2_d_a_d_a = b * b2plusa2_inv2 * (-1) * 2 * a
        d2_atan2_d_a_d_b = (+1) * b2plusa2_inv + b * (-1) * b2plusa2_inv2 * 2 * b
        d2_atan2_d_b_d_a = (-1) * b2plusa2_inv + a * (-1) * b2plusa2_inv2 * 2 * a
        d2_atan2_d_b_d_b = - d2_atan2_d_a_d_a
        
        for i in range(4):
            for j in range(4):
                H[i, j] = (
                        d2_atan2_d_a_d_a * da_dUni[j] * da_dUni[i] + 
                        d2_atan2_d_a_d_b * da_dUni[j] * db_dUni[i] +
                        d2_atan2_d_b_d_a * db_dUni[j] * da_dUni[i] + 
                        d2_atan2_d_b_d_b * db_dUni[j] * db_dUni[i] )
                
        return G, H
    

    
    def G_and_H_Rot_3D(self, Un, Urp, PhiRp, d0):
        
        def Rphi_3D_x( phi, ):
            return np.array([
                                [1,            0,               0,],
                                [0,             np.cos(phi),    -np.sin(phi)  ],
                                [0,              np.sin(phi),   +np.cos(phi)  ],])
    
        def Rphi_3D_y( phi, ):
            return np.array([
                                [np.cos(phi),            0,                +np.sin(phi),],
                                [0,                     1,                  0],
                                [-np.sin(phi),           0,                 +np.cos(phi)  ],])
        def Rphi_3D_z( phi, ):
            return np.array([
                                [ np.cos(phi),  -np.sin(phi),   0  ],
                                [ np.sin(phi),  +np.cos(phi),   0  ],
                                [0,             0,              1,],])
        
        # 3 + 3 + 3
        G = np.zeros( (3, 9) )
        H = np.zeros( (3, 9, 9) )
        
        # d/dUN
        G[:, 0:3, ] = - np.identity(3) 
        # d/dURp
        G[:, 3:6, ] = + np.identity(3) 
        # d
        
        def rot(phi):
            return Rphi_3D_z(phi[2]) @  Rphi_3D_y(phi[1])  @ Rphi_3D_x(phi[0]) 
        
        for i in range(3):
            
            phi = np.copy(PhiRp)
            phi[i] += np.pi/2
            
            G[:, 6 + i] = rot( phi ) @ d0
            
            for j in range(3):
                phi2 = np.copy(phi)
                phi2[j] += np.pi/2
                
                H[:, 6 + i, 6 + j ] = rot( phi2 ) @ d0
        
        return G, H
    
    def applyConstraint(self, Un1, PExt, V, increment):
        
        nConstraints = self.nConstraints
        
        nDim = self.nDim
        
        nU = self.nDof - nConstraints
        nNodes = len ( self.slaveNodes )
        
        UNodes = Un1[: nNodes * nDim]
        
        
        idcsURp = [j for j in range ( nNodes * nDim,  nNodes * nDim + nDim  )]
        idcsPhiRp = [j for j in range ( nNodes * nDim + nDim, nNodes * nDim + nDim + nDim  )]
        

        URp =    Un1[idcsURp]
        PhiRp =  Un1[idcsPhiRp]
        
        Lambdas = Un1[nU:] 
        
        K = V.reshape( self.nDof, self.nDof)
        
        KUU = K[0:nU, 0:nU]
        
        
        for i in range( nNodes ):
        
            d0 = self.d0s[i]
            
            idcsUNode = [ j for j in range ( i * nDim , i * nDim + nDim ) ]
            
            idxL =  (i * nDim )
            
            Lambda = Lambdas[ idxL : idxL + nDim  ]
            
            Un = UNodes [idcsUNode]
                
            G, H = self.G_and_H_Rot_3D(Un, URp, PhiRp, d0)
                
            idcsU = idcsUNode + idcsURp + idcsPhiRp
                                
            KUL = K[0:nU, nU + idxL : nU + idxL + nDim ]
            KLU = K[nU + idxL  : nU + idxL + nDim , 0:nU ]
                
            PExt[idcsU] -= Lambda.T @ G
                
            for k in range(nDim) :
                for i in range(6, 9):
                    for j in range(6, 9):
                        KUU [ idcsU[i], idcsU[j] ]  += Lambda[k] * H[k,i,j]
                    
            KUL [ idcsU, : ] += G.T
            KLU [ :, idcsU ] += G

            # print(K != 0.0)
                    

