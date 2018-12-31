#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 09:18:25 2018

@author: matthias
"""

from collections import OrderedDict
from fe.config.phenomena import getFieldSize
import numpy as np

class VIJSystemMatrix(np.ndarray):
    """This class represents a constructed VIJ triple for sparse matrices in coordinate format,
    and provides direct access to the DOF indices of elements and constraints within the VIJ triple
    """
    
    def __new__(cls, nDof, I, J, entitiesInVIJ ):

        obj = np.zeros_like(I, dtype=float).view(cls)
        
        obj.nDof = nDof
        obj.I = I
        obj.J = J
        obj.entitiesInVIJ = entitiesInVIJ
        
        return obj
    
    def __getitem__(self, key ):
        try:
            idxInVIJ  = self.entitiesInVIJ [ key ]
            return super().__getitem__(  slice (idxInVIJ ,  idxInVIJ + key.nDof ** 2  )  )
        except :
            return super().__getitem__(  key  )
        
class DofVector (np.ndarray) :
    
    def __new__(cls, nDof, entitiesInDofVector ):
        obj = np.zeros(nDof, dtype=float).view(cls)
        obj.entitiesInDofVector = entitiesInDofVector
    
        return obj
    
    def __getitem__(self, key ):
        try:
            return super().__getitem__( self.entitiesInDofVector[key] )
        except :
            return super().__getitem__(key)
    
    def __setitem__(self, key, value):
        try:
            return super().__setitem__( self.entitiesInDofVector[key], value)
        except :
            return super().__setitem__(key, value)    

class DofManager:
    """ The DofManager 
        - analyzes the domain (nodes and constraints), 
        - collects information about the necessary structure of the degrees of freedom 
        - handles the active fields on each node 
        - counts the accumulated number of associated elements on each dof (for the Abaqus like convergence test)
        - supplies the framework with DofVectors and VIJSystemMatrices
        """
    
    def __init__(self, modelInfo):
        
        self.modelInfo = modelInfo
        
        self.activateFieldsOnNodes()
        
        (self.nDof, 
         self.IndicesOfFieldsInDofVector) = self.initializeDofVectorStructure()

        (self.nAccumulatedNodalFluxes, 
         self.nAccumulatedNodalFluxesFieldwise) = self.countNodalFluxes()
        
        (self.accumulatedElementNDof, 
         self.accumulatedConstraintNDof, 
         self.sizeVIJ, 
         self.largestNumberOfElNDof) = self.countAccumulatedEntityDofs()
        
        self.entitiesInDofVector = self.locateEntitiesInDofVector()
        
        (self.I, 
         self.J, 
         self.entitiesInVIJ) = self.initializeVIJPattern()
      
        
    def activateFieldsOnNodes(self, ):
        """ activate all fields on nodes, which are required in the analysis"""
        modelInfo = self.modelInfo
        
        for element in modelInfo['elements'].values():
            for node, nodeFields in zip ( element.nodes, element.fields ):
                node.fields.update( [ (f, True) for f in nodeFields ]  )
                
        for constraint in modelInfo['constraints'].values():
            for node, nodeFields in zip(constraint.nodes, constraint.fieldsOfNodes):
                node.fields.update( [ (f, True) for f in nodeFields]  )
                
    def initializeDofVectorStructure(self):
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
        
        
        IndicesOfFieldsInDofVector = OrderedDict()
        fieldIdxBase = 0
        
        for node in nodes.values():
                #delete all fields of a node, which are not active
            for field, enabled in list(node.fields.items()):
                if not enabled:
                    del node.fields[field]
                    
            for field in node.fields.keys():
                fieldSize = getFieldSize(field, domainSize)
                node.fields[field] = [i + fieldIdxBase for i in range(fieldSize)]
                indexList = IndicesOfFieldsInDofVector.setdefault(field, []) 
                indexList += (node.fields[field])
                fieldIdxBase += fieldSize       
            
        for constraint in constraints.values():
            # some constraints may need additional Degrees of freedom (e.g. lagrangian multipliers)
            # we create them here, and assign them directly to the constraints 
            # (In contrast to true field indices, which are not directly 
            # assigned to elements/constraints but to the nodes)
            nNeededDofs = constraint.getNumberOfAdditionalNeededDofs()
            indicesOfConstraintAdditionalDofs = [i + fieldIdxBase for i in range(nNeededDofs)  ]
            #TODO: Store here in dofmanager
            constraint.assignAdditionalGlobalDofIndices ( indicesOfConstraintAdditionalDofs )
            fieldIdxBase += nNeededDofs
            
        for field, indexList in IndicesOfFieldsInDofVector.items():
            IndicesOfFieldsInDofVector[field] = np.array(indexList)
            
        nDof = fieldIdxBase
        
        return nDof, IndicesOfFieldsInDofVector
                
    def countNodalFluxes(self):
        """
        for the VIJ (COO) system matrix and the Abaqus like convergence test, 
        the number of dofs 'entity-wise' is needed:
        = Σ_(elements+constraints) Σ_nodes ( nDof (field) )"""
                
        IndicesOfFieldsInDofVector = self.IndicesOfFieldsInDofVector
        
        elements    = self.modelInfo['elements']
        constraints = self.modelInfo['constraints']
        
        accumulatedNumberOfFieldFluxes = {}
        accumulatedNodalFluxesTotal = 0
        
        for field in IndicesOfFieldsInDofVector.keys():
             accumulatedNumberOfFieldFluxes[field] =  np.sum(
                     [ len(node.fields[field]) for el in elements.values() 
                                                 for node in (el.nodes) if field in node.fields]
                     + 
                     [ len(node.fields[field]) for constraint in constraints.values() 
                                                 for node in (constraint.nodes) if field in node.fields]) 
             
             accumulatedNodalFluxesTotal += accumulatedNumberOfFieldFluxes[field]
             
        return accumulatedNodalFluxesTotal, accumulatedNumberOfFieldFluxes, 
    
    def countAccumulatedEntityDofs(self):
        """ generates some auxiliary information, 
        which may be required by some modules of EdelweissFE"""

        elements    = self.modelInfo['elements']
        constraints = self.modelInfo['constraints']
        
        sizeVIJ = 0
        accumulatedElementNDof = 0
        largestNumberOfElNDof = 0
        
        for el in elements.values():
            accumulatedElementNDof += el.nDof
            sizeVIJ += (el.nDof**2)
            
            largestNumberOfElNDof =  max ( el.nDof, largestNumberOfElNDof )
        
        nNDofAccumulatedConstraints = 0
        for constraint in constraints.values():
            nNDofAccumulatedConstraints += constraint.nDof
            
            sizeVIJ += constraint.nDof ** 2
            
        return accumulatedElementNDof, nNDofAccumulatedConstraints, sizeVIJ, largestNumberOfElNDof
    
    def locateEntitiesInDofVector(self, ):
        """ Creates a dictionary containing the location (indices) of each entity (elements, constraints)
        within the DofVector structure"""
        
        elements    = self.modelInfo['elements']
        constraints = self.modelInfo['constraints']
        entitiesInDofVector = {}
        
        for el in elements.values():
            destList = np.asarray([i for iNode, node in enumerate(el.nodes) # for each node of the element..
                                        for nodeField in el.fields[iNode]  # for each field of this node
                                            for i in node.fields[nodeField]])  # the index in the global system
    
            entitiesInDofVector[el] = destList[ el.dofIndicesPermutation  ]
            
        for constraint in constraints.values():
            destList = constraint.globalDofIndices
            
            entitiesInDofVector[constraint] = destList
        
        return entitiesInDofVector
                
    def initializeVIJPattern(self, ):
        """ Generate the IJ pattern for VIJ (COO) system matrices"""
        
        elements    = self.modelInfo['elements']
        constraints = self.modelInfo['constraints']
        
        entitiesInVIJ = {}
        entitiesInDofVector = self.entitiesInDofVector
        
        sizeVIJ = self.sizeVIJ
        
        I = np.zeros(sizeVIJ, dtype=np.int)
        J = np.zeros(sizeVIJ, dtype=np.int)
        idxInVIJ = 0
        
        for el in elements.values():
            destList = entitiesInDofVector[el]
    
            entitiesInVIJ[el] = idxInVIJ
                                  
            # looks like black magic, but it's an efficient way to generate all indices of Ke in K:
            elDofLocations = np.tile(destList, (destList.shape[0], 1) )
            I[idxInVIJ : idxInVIJ+el.nDof**2] = elDofLocations.ravel()
            J[idxInVIJ : idxInVIJ+el.nDof**2] = elDofLocations.ravel('F')
            idxInVIJ += el.nDof**2
            
        for constraint in constraints.values():
            destList = entitiesInDofVector[constraint]
            
            entitiesInVIJ[constraint] = idxInVIJ
            
            constraintDofLocations = np.tile( destList, (destList.shape[0], 1) )
            I[idxInVIJ : idxInVIJ + constraint.nDof**2] = constraintDofLocations.ravel()
            J[idxInVIJ : idxInVIJ + constraint.nDof**2] = constraintDofLocations.ravel('F')
            idxInVIJ += constraint.nDof**2
        
        
        return I, J, entitiesInVIJ
    
    def constructVIJSystemMatrix(self, ):
        """ Construct a VIJ (COO) Sparse System matrix object, which also has knowledge about
        the location of each entity"""
        
        nDof            = self.nDof
        I               = self.I
        J               = self.J
        entitiesInVIJ   = self.entitiesInVIJ
        
        return VIJSystemMatrix(nDof, I, J, entitiesInVIJ)
    
    def constructDofVector(self, ):
        """ Construct a vector with size=nDof and which has knowledge about 
        the location of each entity"""
        
        nDof = self.nDof
        entitiesInDofVector = self.entitiesInDofVector
        
        return DofVector(nDof, entitiesInDofVector)