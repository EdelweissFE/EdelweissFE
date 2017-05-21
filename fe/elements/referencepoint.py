#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 13:42:05 2017

@author: matthias
"""

from fe.config.phenomena import phenomena
from collections import OrderedDict

class ReferencePoint:
    """ Base class for a finite element refernce Point,
    currently nothing more than a dictioniary for label and coordinates"""
    
    def __init__(self, label, coordinates,):
        
        self.label = label
        self.coordinates = coordinates
        self.fields = OrderedDict.fromkeys(phenomena, False)
    
    def setFields(self, *fields):
        self.fields.update(fields)
    