#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 19:20:25 2017

@author: matthias
"""

import importlib

solverLibrary = {'NIST' : 'nonlinearimplicitstatic',
                 'NISTParallel' : 'nonlinearimplicitstaticparallel',
                 'NISTPArcLength' : 'nonlinearimplicitstaticparallelarclength',
                 }


def getSolverByName(name):
    solver = importlib.import_module("fe.solvers.{:}".format( solverLibrary[name] ))    
    return getattr ( solver, name ) 

