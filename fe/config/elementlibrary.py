#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""
#import importlib

def getElementByName(name):
    
    if name.lower().startswith('uel'):
        from fe.elements.marmotelement.element import MarmotElementWrapper
        return MarmotElementWrapper
