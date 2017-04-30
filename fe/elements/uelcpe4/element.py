import numpy as np
from fe.elements.uelbaseelement.element import BaseElement

class Element(BaseElement):
    fields =                [["displacement"],
                             ["displacement"],
                             ["displacement"],
                             ["displacement"],] # fields identical for each node
    nNodes =                4
    nGaussPt =              4
    nDofPerEl =             8
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation  = np.arange(0, 8, 1)
    ensightType =           "quad4"
    uelIdentification =     407
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.uelIdentification)
        
    def setProperties(self, uelProperties, umatName, nStateVarsUmat):
        super().setProperties(uelProperties, umatName, nStateVarsUmat)
                
