#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 11:02:12 2017

@author: matthias
"""

class WrongDomain(Exception):
    pass

class StepFailed(Exception):
    pass

class CutbackRequest(Exception):
    def __init__(self, message, cutbackSize):
        super().__init__(message)
        self.cutbackSize = float(cutbackSize)
        
class ReachedMaxIterations(Exception):
    pass
class ReachedMaxIncrements(Exception):
    pass
class ReachedMinIncrementSize(Exception):
    pass
