#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""

# plane strain continuum elements
from  fe.elements.uelcpe4.element import Element as uelCPE4
from  fe.elements.uelcpe4r.element import Element as uelCPE4R

# plane stress continuum elements
# 4 node
from  fe.elements.uelcps4.element import Element as uelCPS4
from  fe.elements.uelcps4nonlocal.element import Element as uelCPS4NonLocal
## 8 node
from  fe.elements.uelcps8r.element import Element as uelCPS8R
from  fe.elements.uelcps8rnonlocal.element import Element as uelCPS8RNonLocal
from  fe.elements.uelcps8nonlocal.element import Element as uelCPS8NonLocal
from  fe.elements.uelcpe8rnonlocal.element import Element as uelCPE8RNonLocal

elementlibrary = {  
                    "uelCPE4" :         uelCPE4,
                    "uelCPE4R" :        uelCPE4R,
                    "uelCPS4" :         uelCPS4,
                    "uelCPS4NonLocal" : uelCPS4NonLocal,
#                    
                     "uelCPS8R" :        uelCPS8R,
                    "uelCPS8RNonLocal": uelCPS8RNonLocal,
                    "uelCPS8NonLocal": uelCPS8NonLocal,
                    # "uelCPE8RNonLocal": uelCPE8RNonLocal,
                    "uelCPS8RNonLocal": uelCPS8RNonLocal,
                    "uelCPE8RNonLocal": uelCPE8RNonLocal,
                    }
