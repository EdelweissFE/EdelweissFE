import numpy as np
from fe.elements.uelbaseelement.element import BaseElement

class Element(BaseElement):
    fields =                [["displacement", "nonlocal damage"], # node 1
                             ["displacement", "nonlocal damage"], # node 2
                             ["displacement", "nonlocal damage"], # node 3
                             ["displacement", "nonlocal damage"], # node 4
                             ["displacement", "nonlocal damage"], # node 5
                             ["displacement", "nonlocal damage"], # node 6
                             ["displacement", "nonlocal damage"], # node 7
                             ["displacement", "nonlocal damage"],]# node 8
    nNodes =                8
    nGaussPt =              9
    nDofPerEl =             24
    sizeKe =                24 * 24
    dofIndicesPermutation = np.array([0,1,3,4,6,7,9,10,12,13,15,16,18,19,21,22] + [2,5,8,11,14,17,20,23], dtype=int)   
    ensightType =           "quad8"
    uelIdentification =     817
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.uelIdentification)
