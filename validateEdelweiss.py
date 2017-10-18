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
    parser.add_argument('--create', dest='create', action='store_true', help='create refernce solutions') 
    args=parser.parse_args()

    testfile =  'test.inp' 
    referenceSolutionFile = 'U.ref'
    
    testfilesDir = os.path.join ( os.getcwd() , 'testfiles' )
    os.chdir(testfilesDir)
    testsDirs = next(os.walk('.'))[1]
    
    for directory in testsDirs:
        
        os.chdir( os.path.join ( testfilesDir, directory ))
        
        if not os.path.exists(testfile):
            continue
        
        inputFile = parseInputFile(testfile)
        print('Test {:30}'.format(directory ), end='\r')
        
        try:
        
            success, U, P, fieldOutputController = finitElementSimulation(inputFile, verbose = False, suppressPlots=True)
            
            if not args.create:
                
                UReference = np.loadtxt(referenceSolutionFile)
                residual = U - UReference
                
                if success and ( np.max ( np.abs ( residual ) ) ) < 1e-6:
                    passed = True
                else:
                    passed = False
                
                if passed:
                    print('Test {:30} PASSED'.format(directory))
                else:
                    print('Test {:30} FAILED'.format(directory))
                os.chdir('..')
                
            else:
                print('')
                np.savetxt(referenceSolutionFile, U)

        except Exception as e:
            print("Test {:30} FAILED: ".format(directory) + str(e) )
            continue