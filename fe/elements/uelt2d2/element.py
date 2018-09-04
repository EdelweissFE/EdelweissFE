#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 14 19:22:23 2018

@author: matthias
"""

from fe.elements.uelbaseelement.element import BaseElement


class Element(BaseElement):
    fields =                [["displacement", ], # node 1
                             ["displacement", ], # node 2
                             ]
    
    nNodes =                    2
    nGaussPt =                  2
    nDofPerEl =                 4
    sizeKe =                    nDofPerEl * nDofPerEl
    dofIndicesPermutation =     slice(0, 4)
    ensightType =               "bar2"
    uelIdentification =         "UelT2D2"
    nStateVarsGaussPtSpecific = 12
    
    def __init__(self, nodes, elNumber):
        super().__init__(nodes, elNumber, self.uelIdentification)
