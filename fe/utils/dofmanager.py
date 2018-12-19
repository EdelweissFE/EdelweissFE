#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 09:18:25 2018

@author: matthias
"""

from collections import OrderedDict
from fe.config.phenomena import getFieldSize
import numpy as np

class VIJDatabase:
    """This class represents a constructed VIJ triple for sparse matrices in coordinate format,
    and provides direct access to the DOF indices of elements and constraints within the VIJ triple
    """
    
    def __init__(self, elements, constraints):
        
        entitiesInVIJ = {}# self.entitiesInVIJ
        
        self.sizeVIJ = 0
        sizeVIJ = self.sizeVIJ
        
        for el in elements.values():
            sizeVIJ += (el.nDofPerEl**2)
        
        for constraint in constraints.values():
            sizeVIJ += constraint.sizeStiffness
        
        V = np.zeros(sizeVIJ)
        I = np.zeros_like(V, dtype=np.int)
        J = np.zeros_like(V, dtype=np.int)
        idxInVIJ = 0
        
        for el in elements.values():
            destList = np.asarray([i for iNode, node in enumerate(el.nodes) # for each node of the element..
                                        for nodeField in el.fields[iNode]  # for each field of this node
                                            for i in node.fields[nodeField]])  # the index in the global system
    
            entitiesInVIJ[el] = idxInVIJ
                                  
            # looks like black magic, but it's an efficient way to generate all indices of Ke in K:
            elDofLocations = np.tile(destList[ el.dofIndicesPermutation  ], (destList.shape[0], 1) )
            I[idxInVIJ : idxInVIJ+el.nDofPerEl**2] = elDofLocations.ravel()
            J[idxInVIJ : idxInVIJ+el.nDofPerEl**2] = elDofLocations.ravel('F')
            idxInVIJ += el.nDofPerEl**2
            
        for constraint in constraints.values():
            destList = constraint.globalDofIndices
            
            entitiesInVIJ[constraint] = idxInVIJ
            
            constraintDofLocations = np.tile( destList, (destList.shape[0], 1) )
            I[idxInVIJ : idxInVIJ + constraint.sizeStiffness] = constraintDofLocations.ravel()
            J[idxInVIJ : idxInVIJ + constraint.sizeStiffness] = constraintDofLocations.ravel('F')
            idxInVIJ += constraint.sizeStiffness
        
        self.V = V
        self.I = I
        self.J = J
        self.entitiesInVIJ = entitiesInVIJ

class DofManager:
    """ The DofManager 
        - analyzes the domain (nodes and constraints), 
        - collects information about the necessary structure of the degrees of freedom 
        - handles the active fields on each node 
        - counts the accumulated number of associated elements on each dof (for the Abaqus like convergence test)
        """
    
    def __init__(self, modelInfo):
        
        self.modelInfo = modelInfo
        
        self.createFieldsOnNodes()
        
        self.numberOfDofs, self.fieldIndices = self.assignFieldDofIndices()
        
        self.numberOfAccumulatedFieldConnections = self.countAccumulatedFieldConnections()
      
        
    def createFieldsOnNodes(self, ):
        
        modelInfo = self.modelInfo
        
        for element in modelInfo['elements'].values():
            for node, nodeFields in zip ( element.nodes, element.fields ):
                node.fields.update( [ (f, True) for f in nodeFields ]  )
                
        for constraint in modelInfo['constraints'].values():
            for node, nodeFields in zip(constraint.nodes, constraint.fieldsOfNodes):
                node.fields.update( [ (f, True) for f in nodeFields]  )
                
    def assignFieldDofIndices(self):
        """ Loop over all nodes to generate the global field-dof indices.
        output is a tuple of:
            - number of total DOFS
            - orderedDict( (mechanical, indices), 
                           (nonlocalDamage, indices)
                           (thermal, indices)
                           ...)."""
        
        nodes = self.modelInfo['nodes']
        domainSize = self.modelInfo['domainSize']
        constraints = self.modelInfo['constraints']
        
        
        fieldIndices = OrderedDict()
        fieldIdxBase = 0
        
        for node in nodes.values():
                #delete all fields of a node, which are not active
            for field, enabled in list(node.fields.items()):
                if not enabled:
                    del node.fields[field]
                    
            for field in node.fields.keys():
                fieldSize = getFieldSize(field, domainSize)
                node.fields[field] = [i + fieldIdxBase for i in range(fieldSize)]
                indexList = fieldIndices.setdefault(field, []) 
                indexList += (node.fields[field])
                fieldIdxBase += fieldSize             
            
        for constraint in constraints.values():
            # some constraints may need additional Degrees of freedom (e.g. lagrangian multipliers)
            # we create them here, and assign them directly to the constraints 
            # (In contrast to true field indices, which are not directly 
            # assigned to elements/constraints but to the nodes)
            nNeededDofs = constraint.getNumberOfAdditionalNeededDofs()
            indicesOfConstraintAdditionalDofs = [i + fieldIdxBase for i in range(nNeededDofs)  ]
            constraint.assignAdditionalGlobalDofIndices ( indicesOfConstraintAdditionalDofs )
            fieldIdxBase += nNeededDofs
            
        for field, indexList in fieldIndices.items():
            fieldIndices[field] = np.array(indexList)
            
        numberOfDofs = fieldIdxBase
        
        return numberOfDofs, fieldIndices
                
    def countAccumulatedFieldConnections(self):
        """
        for the Abaqus like convergence test, the number of dofs 'element-wise' is needed:
        = Σ_(elements+constraints) Σ_nodes ( nDof (field) )"""
                
        fieldIndices = self.fieldIndices
        
        elements = self.modelInfo['elements']
        constraints = self.modelInfo['constraints']
        
        accumulatedNumberOfFieldConnections = {}
        
        for field in fieldIndices.keys():
             accumulatedNumberOfFieldConnections[field] =  np.sum(
                     [ len(node.fields[field]) for el in elements.values() 
                                                 for node in (el.nodes) if field in node.fields]
                     + 
                     [ len(node.fields[field]) for constraint in constraints.values() 
                                                 for node in (constraint.nodes) if field in node.fields]) 
        return accumulatedNumberOfFieldConnections
                
    def constructVIJDatabase(self, ):
        """ Initializes the V vector and generates I, J entries for each element,
        based on i) its (global) nodes ii) its dofLayout. Furthermore, 
        a dictionary with the mapping of each element to its index in VIJ 
        is created.
        The same is done for each constraint
        -> is called by __init__() """
        
        elements    = self.modelInfo['elements']
        constraints = self.modelInfo['constraints']
        
        return VIJDatabase(elements, constraints)