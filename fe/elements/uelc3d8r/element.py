import numpy as np
from fe.elements.uelbaseelement.element import BaseElement

class Element(BaseElement):
    fields =                [["displacement", ], # node 1
                             ["displacement", ], # node 2
                             ["displacement", ], # node 3
                             ["displacement", ], # node 4
                             ["displacement", ], # node 5
                             ["displacement", ], # node 6
                             ["displacement", ], # node 7
                             ["displacement", ], # node 8
                             ]
    
    nNodes =                8
    nGaussPt =              8
    nDofPerEl =             24
    sizeKe =                24 * 24
    dofIndicesPermutation = slice(0, 24)#np.arange(0, 24, 1, dtype=np.int)
    ensightType =           "hexa8"
    uelIdentification =     "UelC3D8R"
    nStateVarsGaussPtSpecific =     12
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.nStateVarsGaussPtSpecific, self.uelIdentification)
