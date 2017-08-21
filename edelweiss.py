#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""

import argparse  
from fe.fecore import finitElementSimulation
from fe.parameteridentification.paramIdentification import parameterIdentification
from fe.utils.inputfileparser import parseInputFile,printKeywords
from fe.utils.printdocumentation import printDocumentation

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='Batch computation for Edelweiss finite element jobs')
    
    parser.add_argument('file', type=str,  nargs='*', ) # multiple input files possible
    parser.add_argument('--quiet', dest='verbose', action='store_false', help='suppress output') 
    parser.add_argument('--keywords', dest='kw', action='store_true', help='print keywords') 
    parser.add_argument('--doc=module', dest='doc', help='print keywords') 
    parser.add_argument('--parameterIdentification', dest='parID', help='run parameter identification',  nargs='*',)
    args=parser.parse_args()
    
    fileList=args.file

    if args.kw:
        printKeywords()
        exit(0)
    
    if args.doc:
        printDocumentation(args.doc)
        exit(0)
    
    inputFiles = []
    #1 ) parse all files
    try:
        for file in fileList:    
          #  print(parseInputFile(file))
            inputFiles.append(parseInputFile(file))

    except (KeyError, ValueError) as e:
        print(str(e))
        exit(1)
    
    if args.parID:
        for inputFile in inputFiles:
            parameterIdentification(inputFile)
            
#    #2 ) all computations and imports
    else:
        for inputFile in inputFiles:
            finitElementSimulation(inputFile, verbose = args.verbose)
