import numpy as np
from fe.elements.uelbaseelement.element import BaseElement

class Element(BaseElement):
    fields =                [["displacement", ], # node 1
                             ["displacement", ], # node 2
                             ["displacement", ], # node 3
                             ["displacement", ],]# node 4
    nNodes =                4
    nGaussPt =              4
    nDofPerEl =             8
    sizeKe =                8 * 8
    dofIndicesPermutation = np.array([0,1,2,3,4,5,6,7])
    ensightType =           "quad4"
    uelIdentification =     402
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.nGaussPt, self.uelIdentification)
        
    def setProperties(self, uelProperties, umatName, nStateVarsUmat):
        super().setProperties(uelProperties, umatName, nStateVarsUmat)
                
