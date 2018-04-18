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
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],
                             ["displacement", "nonlocal damage"],] 
    
    nNodes =                20
    nGaussPt =              8
    nDofPerEl =             nNodes * 4
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation = np.asarray( [  k for i in range(nNodes) for k in [i*4, (i*4+1), (i*4+2)]   ] + [ (i*4 + 3 ) for i in range (nNodes)] , dtype=int) 
    ensightType =           "hexa20"
    uelIdentification =     "UelC3D20RNonLocal"
    nStateVarsGaussPtSpecific =     12
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.nStateVarsGaussPtSpecific, self.uelIdentification)
