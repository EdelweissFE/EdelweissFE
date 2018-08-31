#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 08:50:46 2017

@author: matthias
"""

import numpy as np

def sigPrinc(x):
    return np.linalg.eig([[x[0], x[3], x[4]],
                          [x[3], x[1], x[5]],
                          [x[4], x[5], x[2]]] )[0]

mathModules = {'sigPrinc':sigPrinc,
               'mean'   : np.mean,
               'max'    : np.max,
               'min'    : np.min,
               'linalg' : np.linalg}

def createModelAccessibleFunction(expression, modelInfo, fieldOutputs):
    """ Create a function from a string expression, which can access the complete model and fieldOutput"""
    scope = {**modelInfo, **locals()}
    f = eval("lambda:" + expression, scope)
    return f
    
def evalModelAccessibleExpression(expression, modelInfo, fieldOutputs):
    """ Evalualate a string expression, which can access the complete model and fieldOutput"""
    return createModelAccessibleFunction(expression, modelInfo, fieldOutputs)()

def createMathExpression(expression, symbol='x'):
    return lambda x : eval ( expression, globals(), { symbol : x , **mathModules } )
