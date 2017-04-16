#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""
import importlib

def getElementByName(name):
    module = importlib.import_module("fe.elements."+name.lower()+".element")    
    return module.Element