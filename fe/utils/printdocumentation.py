#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 19:38:14 2017

@author: matthias
"""
import importlib
import textwrap

def printDocumentation(module):
    
    try:
        mod = importlib.import_module ('fe.'+module )
        print(mod.__doc__)
        wrapper = textwrap.TextWrapper(width=80,replace_whitespace=False)
        for key, val in mod.documentation.items():
            wrapper.initial_indent = "    {:26}".format(key)
            wrapper.subsequent_indent = ' '*len(wrapper.initial_indent)
            print(wrapper.fill(val))
            
        
    except ModuleNotFoundError as e:
        print(e)
    except AttributeError as e:
        print('documentation for module {:} not found'.format(module))
    