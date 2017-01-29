#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 19:53:45 2017

@author: matthias
"""

from fe.config.phenomena import phenomena
from collections import OrderedDict

class Node:
    """ Base class for a finite element node,
    currently nothing more than a dictioniary for label and coordinates"""
    
    def __init__(self, label, coordinates,):
        
        self.label = label
        self.coordinates = coordinates
        self.fields = OrderedDict.fromkeys(phenomena, False)
    
    def setFields(self, *fields):
        self.fields.update(fields)
    