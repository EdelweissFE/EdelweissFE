#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:37:28 2017

@author: matthias
"""

import importlib

def getGeneratorByName(name):
    module = importlib.import_module("fe.generators."+name.lower())    
    return module.generateModelData