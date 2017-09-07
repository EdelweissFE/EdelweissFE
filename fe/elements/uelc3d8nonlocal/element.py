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
                             ["displacement", "nonlocal damage"],] 
    
    nNodes =                8
    nGaussPt =              8
    nDofPerEl =             32
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation = np.asarray([0, 1, 2, 4, 5, 6, 8, 9, 
                                        10, 12, 13, 14, 16, 17, 18, 20,
                                        21, 22, 24, 25, 26, 28, 29, 30, 
                                        3, 7,11, 15, 19, 23, 27, 31], dtype=int)
    ensightType =           "hexa8"
    uelIdentification =     813
    nStateVarsGaussPtAdditional =     12
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.nStateVarsGaussPtAdditional, self.uelIdentification)