#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 19:20:25 2017

@author: matthias
"""

from fe.solvers.nonlinearimplicitstatic import NIST
from fe.solvers.nonlinearimplicitstaticparallel import NISTParallel

solverLibrary = {'NIST' : NIST,
                 'NISTParallel' : NISTParallel,}