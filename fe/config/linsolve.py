#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 10:27:25 2018

@author: matthias
"""

def getDefaultLinSolver():
    try:
        from fe.linsolve.pardiso.pardiso import pardisoSolve
        return pardisoSolve
    except:
        from scipy.sparse.linalg import spsolve 
        return lambda A, b : spsolve(A, b, use_umfpack=False)

def getLinSolverByName(linsolverName):

    if linsolverName.lower() == 'superlu':
        from scipy.sparse.linalg import spsolve 
        return lambda A, b : spsolve(A, b, use_umfpack=False)
    elif linsolverName.lower() == 'umfpack':
        from scipy.sparse.linalg import spsolve 
        return lambda A, b : spsolve(A, b, use_umfpack=False)
    elif linsolverName.lower() == 'pardiso':
        from fe.linsolve.pardiso.pardiso import pardisoSolve
        return pardisoSolve
    else:
        raise AttributeError('invalid linear solver {:} requested'.format(linsolverName))