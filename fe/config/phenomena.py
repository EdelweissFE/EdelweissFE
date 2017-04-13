#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  10 19:10:42 2017

@author: matthias
"""
from collections import OrderedDict

                                        #field              field type
phenomena = OrderedDict([               ("mechanical" ,     "vector"),
                                        ("thermal",         "scalar"),
                                        ("nonlocal damage", "scalar"),])
                                        
                                        
                                        #field              #tolerance
flowCorrectionTolerance = {             'mechanical' :      1e-5,
                                        'nonlocal damage':  1e-5,
                                        }
                                                                   
effortResidualTolerance = {             'mechanical' :      1e-7,
                                        'nonlocal damage':  1e-7,
                                        }

effortResidualToleranceAlternative = {  'mechanical' :       1e-3,
                                        'nonlocal damage':   1e-6,
                                        }
                                        
                                        #domain             #dimensions
domainMapping = {                       "1d":               1,
                                        "2d":               2,
                                        "3d":               3,
                                        "axisymmetric":     2,
                                        }
    
def getFieldSize(field, domainSize):
    fType = phenomena[field]
    if fType == "scalar":
        return 1
    if fType == "vector":
        return domainSize

    
