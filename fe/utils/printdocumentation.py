#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 19:38:14 2017

@author: matthias
"""
import importlib

def printDocumentation(module):
    
    try:
        mod = importlib.import_module ('fe.'+module )
        mod.printDocumentation()
        
    except ModuleNotFoundError as e:
        print(e)
    except AttributeError as e:
        print('documentation for module {:} not found'.format(module))
    