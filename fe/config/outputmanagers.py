#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 09:27:14 2017

@author: matthias
"""
import importlib

def getOutputManagerByName(name):
    module = importlib.import_module("fe.outputmanagers."+name.lower())    
    return module.OutputManager
