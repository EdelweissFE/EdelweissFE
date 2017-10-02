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
                             ["displacement", ], # node 9
                             ["displacement", ], # node 10
                             ["displacement", ], # node 11
                             ["displacement", ], # node 12
                             ["displacement", ], # node 13
                             ["displacement", ], # node 14
                             ["displacement", ], # node 15
                             ["displacement", ], # node 16
                             ["displacement", ], # node 17
                             ["displacement", ], # node 18
                             ["displacement", ], # node 19
                             ["displacement", ], # node 20
                             ]
    
    nNodes =                len(fields)
    nGaussPt =              27
    nDofPerEl =             nNodes * 3
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation = slice(0, nDofPerEl)
    ensightType =           "hexa20"
    uelIdentification =     "UelC3D20"
    nStateVarsGaussPtAdditional =     12
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.nStateVarsGaussPtAdditional, self.uelIdentification)
