import numpy as np
from fe.elements.uelbaseelement.element import BaseElement

class Element(BaseElement):
    fields =                [["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],] # fields identical for each node
    
    nNodes =                8
    nGaussPt =              4
    nDofPerEl =             24
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation = np.array([0,1,3,4,6,7,9,10,12,13,15,16,18,19,21,22] + [2,5,8,11,14,17,20,23], dtype=int)  
    ensightType =           "quad8"
    uelIdentification =     "UelCPE8RNonLocal"
    nStateVarsGaussPtSpecific =     12
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.nStateVarsGaussPtSpecific, self.uelIdentification)