#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 13:05:53 2017

@author: matthias
"""

import importlib

def getStepActionGeneratorByName(name):
    module = importlib.import_module("fe.stepactions."+name.lower())    
    return module.generateAction