#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""

import argparse  
from fe.fecore import finitElementSimulation
from fe.utils.inputfileparser import parseInputFile

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='Batch computation for Edelweiss finite element jobs')
    parser.add_argument('file', type=str,  nargs='+', ) # multiple input files possible
    parser.add_argument('--quiet', dest='verbose', action='store_false') 
    args=parser.parse_args()
    
    fileList=args.file
    inputFiles = []
    #1 ) parse all files
    try:
        for file in fileList:     
                inputFiles.append(parseInputFile(file))
    except (KeyError, ValueError) as e:
        print(e)
        exit(1)
        
    #2 ) all computations and imports
    for inputFile in inputFiles:
        finitElementSimulation(inputFile, verbose = args.verbose)
