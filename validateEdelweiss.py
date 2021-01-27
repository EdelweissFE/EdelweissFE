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
import argparse  

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='validation script for FE analyses')
    parser.add_argument('--create', dest='create', action='store_true', help='create refernece solutions') 
    parser.add_argument('--tests', help='comma-separated list (without whitespace inbetween) with names of analyzed test files, ' 
                                         'e.g. MeshPlot,NodeForces, or simply type all. The names are case-sensitive.', type=str, default='all')
    args=parser.parse_args()

    testfile =  'test.inp' 
    referenceSolutionFile = 'U.ref'
     
    tests = [item for item in args.tests.split(',')]

    testfilesDir = os.path.join ( os.getcwd() , 'testfiles' )
    os.chdir(testfilesDir)
    testsDirs = next(os.walk('.'))[1]
    
    if "all" not in tests:
        testsDirs = list(set(testsDirs).intersection(set(tests)))       

    for directory in testsDirs:
        
        os.chdir( os.path.join ( testfilesDir, directory ))
        
        # no test.inp file is found
        if not os.path.exists(testfile):
            continue

        inputFile = parseInputFile(testfile)
        print('Test {:30}'.format(directory ), end='\r')
        
        try:
            success, U, P, fieldOutputController = finitElementSimulation(inputFile, verbose = False, suppressPlots=True)
            
            if not args.create:
                if success:
                    UReference = np.loadtxt(referenceSolutionFile)
                    residual = U - UReference
                    
                    if ( np.max ( np.abs ( residual ) ) ) < 1e-6:
                        print('Test {:30} PASSED'.format(directory))
                    else:
                        print('Test {:30} FAILED'.format(directory))
                    os.chdir('..')
                else:
                    print('Test {:30} FAILED: '.format(directory) + "Test not completed!")
                
            else:
                print('')
                np.savetxt(referenceSolutionFile, U)

        except Exception as e:
            print("Test {:30} FAILED: ".format(directory) + str(e) )
            continue
