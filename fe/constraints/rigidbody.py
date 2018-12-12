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
       
        nConstraints = 0
        rotationalConstraint = True
        if rotationalConstraint:
            if self.nDim == 2 :
                nConstraints += nSlaves * 1 
            elif self.nDim == 3:
                nConstraints += nSlaves * 3
        
        nRp = 1
        
        nRpRotDof = 1 if self.nDim == 2 else 3
        self.nAffectedDofs = self.nDim * (nSlaves + nRp ) + nRpRotDof 

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
    
    def applyConstraint(self, Un1, PExt, V, increment):
        
        nConstraints = self.nConstraints
        
        nDim = self.nDim
        
        if nDim == 2:
            nRot = 1
        elif nDim == 3:
            nRot = 3
        
        nU = self.nDof - nConstraints
        nNodes = len ( self.slaveNodes )
        
        UNodes = Un1[: nNodes * nDim]
        
        
        idcsURp = [j for j in range ( nNodes * nDim,  nNodes * nDim + nDim  )]
        idcsPhiRp = [j for j in range ( nNodes * nDim + nDim, nNodes * nDim + nDim + nRot  )]
        

        URp =    Un1[idcsURp]
        PhiRp =  Un1[idcsPhiRp]
        
        Lambdas = Un1[nU:] 
        
        K = V.reshape( self.nDof, self.nDof)
        
        KUU = K[0:nU, 0:nU]
        
        idcsXY = [0,1]
        idcsXZ = [0,2]
        idcsYZ = [1,2]
        
        for i in range( nNodes ):
        
            d0 = self.d0s[i]
            
            idcsUNode = [ j for j in range ( i * nDim , i * nDim + nDim ) ]
            
            for n in range(nRot):
                #projection in each direction
                if n == 0:
                    coordIndices = idcsXY
                elif n == 1:
                    coordIndices = idcsXZ
                else:
                    coordIndices = idcsYZ
                
                d02 = d0 [ coordIndices ]

                idxL =  (i * nRot ) + n 
                Lambda = Lambdas[ idxL]
                
                Un = UNodes [idcsUNode]
                
                Un2 = Un [ coordIndices ]
                URp2 = URp [coordIndices]
                
                G, H =  self.G_and_H_Rot(Un2[0], Un2[1], URp2[0], URp2[1], PhiRp[n], d02)
                
                idcsU = idcsUNode + idcsURp + idcsPhiRp
                                
                KUL = K[0:nU, nU+ idxL ]
                KLU = K[nU + idxL, 0:nU ]
                
                for j in range(5):
                    PExt[idcsU[j]] -= Lambda * G[j]
                    
                    for k in range(5):
                        KUU [ idcsU[j], idcsU[k] ] += Lambda * H[j,k]
                        
                    KUL [ idcsU[j ]] += G[j]
                    KLU [ idcsU[j] ] += G[j]
                    

