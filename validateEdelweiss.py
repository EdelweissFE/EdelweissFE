#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 21:40:55 2017

@author: matthias
"""

import os

from fe.fecore import finitElementSimulation
from fe.utils.inputfileparser import parseInputFile
import numpy as np

testfile =  'test.inp' 
referenceSolutionFile = 'U.ref'

testfilesDir = os.path.join ( os.getcwd() , 'testfiles' )
os.chdir(testfilesDir)
testsDirs = next(os.walk('.'))[1]

for directory in testsDirs:
    
    os.chdir( os.path.join ( testfilesDir, directory ))
    
    if not os.path.exists(testfile):
        continue
    
    UReference = np.loadtxt(referenceSolutionFile)
    
    inputFile = parseInputFile(testfile)
    print('Test {:30}'.format(directory ), end='\r')
    success, U, P, fieldOutputController = finitElementSimulation(inputFile, verbose = False)
    
    residual = U - UReference
    
    if ( np.max ( np.abs ( residual ) ) ) < 1e-6:
        passed = True
    else:
        passed = False
    
    if passed:
        print('Test {:30} PASSED'.format(directory))
    else:
        print('Test {:30} FAILED'.format(directory))
    os.chdir('..')
