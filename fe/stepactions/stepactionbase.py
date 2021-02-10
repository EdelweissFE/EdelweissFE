#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 20:05:41 2017

@author: matthias
"""

from abc import ABC, abstractmethod

class StepActionBase(ABC):
    """ This is the abstract base class for all step actions managers.
    User defined step actions must implement the abstract methods."""
    
    identification = "StepActionBase"
    
    
    @abstractmethod
    def __init__(self, name, definition, jobInfo, modelInfo, fieldOutputController, journal):
        pass
    
    @abstractmethod
    def updateStepAction(self, definition):
        """is called when an updated definition is present for a new step"""
        pass
    
    @abstractmethod
    def finishStep(self, U, P):
        """is called when a step successfully finished"""
        pass
