#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 21:40:55 2017

@author: matthias
"""

import os

import matplotlib
matplotlib.use('Agg')

from fe.fecore import finitElementSimulation
from fe.utils.inputfileparser import parseInputFile
import numpy as np
import argparse  
from timeit import default_timer as timer
from rich import print
import rich


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
        print('Test {:50}'.format(directory ), end='\r')
        
        try:
       
            tic = timer()
            success, U, P, fieldOutputController = finitElementSimulation(inputFile, verbose = False, suppressPlots=True)
            toc = timer()
            
            if not args.create:
                if success:
                    UReference = np.loadtxt(referenceSolutionFile)
                    residual = U - UReference
                    
                    if ( np.max ( np.abs ( residual ) ) ) < 1e-6:
                        print("Test {:50} [green]PASSED[/] \[{:2.1f}]".format(directory, toc-tic) )
                    else:
                        print('Test {:50} [red]FAILED[/]'.format(directory))
                    os.chdir('..')
                else:
                    print('Test {:50} [red]FAILED[/]: '.format(directory) + "Test not completed!")
                
            else:
                print('')
                np.savetxt(referenceSolutionFile, U)

        except Exception as e:
            print("Test {:50} [red]FAILED[/]: ".format(directory) + str(e) )
            continue
