#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 19:35:09 2017

@author: matthias
"""
from fe.utils.misc import stringDict

def generateNISTSolverOptions(actionDefinitionLines, jobInfo, modelInfo, time, 
                                                               stepActions, 
                                                               U, P):
    oldOptions = stepActions.get('NISTSolverOptions', {})
    newOptions = stringDict([entry for line in actionDefinitionLines for entry in line])
    
    oldOptions.update(newOptions)
    return oldOptions

        