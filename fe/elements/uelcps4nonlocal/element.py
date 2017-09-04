import numpy as np
from fe.elements.uelbaseelement.element import BaseElement

class Element(BaseElement):
    fields =                [["displacement", "nonlocal damage"], # node 1
                             ["displacement", "nonlocal damage"], # node 2
                             ["displacement", "nonlocal damage"], # node 3
                             ["displacement", "nonlocal damage"],]# node 4
    nNodes =                4
    nGaussPt =              4
    nDofPerEl =             12
    sizeKe =                12 * 12
    dofIndicesPermutation = np.array([0,1,3,4,6,7,9,10,2,5,8,11])
    ensightType =           "quad4"
    uelIdentification =     412
    nStateVarsGaussPtAdditional =     12
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.nStateVarsGaussPtAdditional, self.uelIdentification)
