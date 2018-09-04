import numpy as np
from fe.elements.uelbaseelement.element import BaseElement

class Element(BaseElement):
    fields =                [["displacement"],
                             ["displacement"],
                             ["displacement"],
                             ["displacement"],
                             ["displacement"],
                             ["displacement"],
                             ["displacement"],
                             ["displacement"],                              
                             ] # fields identical for each node
    
    nNodes =                8
    nGaussPt =              4
    nDofPerEl =             16
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation  = np.arange(0, 16, 1)
    ensightType =           "quad8"
    uelIdentification =     "UelCPS8R"
    nStateVarsGaussPtSpecific =     12
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.uelIdentification)
