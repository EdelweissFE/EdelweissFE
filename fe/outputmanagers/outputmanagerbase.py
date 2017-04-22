#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:09:35 2017

@author: matthias
"""

from abc import ABC, abstractmethod

class OutputManagerBase(ABC):
    """ This is the abstract base class for all output managers.
    User defined output managers must implement the abstract methods."""
    
    identificatiom = "OutputManager"
    
    @abstractmethod
    def __init__(self, name, definitionLines, jobInfo, modelInfo, journal):
        pass
    
    @abstractmethod
    def initializeStep(self, step, stepActions):
        pass
    
    @abstractmethod
    def finalizeIncrement(self, U, P, increment):
        pass
    
    @abstractmethod
    def finalizeStep(self, U, P):
        pass
    
    @abstractmethod
    def finalizeJob(self, U, P):
        pass