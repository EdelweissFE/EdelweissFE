#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 13:05:53 2017

@author: matthias
"""

from fe.stepactions.dirichlet import generateDirichlet
from fe.stepactions.nodeforces import generateNodeForces
from fe.stepactions.nistsolveroptions import generateNISTSolverOptions

stepActionModules = { 'dirichlet' : generateDirichlet,
                      'nodeForces': generateNodeForces,
                      'NISTSolverOptions': generateNISTSolverOptions,
        }
