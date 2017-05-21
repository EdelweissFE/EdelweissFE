#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:27:14 2017

@author: matthias
"""
import importlib

def getConstraintByName(name):
    module = importlib.import_module("fe.constraints."+name.lower())    
    return module.Constraint
