#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 11:30:59 2017

@author: matthias
"""

import math

class Journal:
    """ This class provides an interface to present messages to the user via console
    output and/or file output.
    Information messages can be sorted by the importance level. 
    Suppressing certain levels of output is possible. """
    
    outputWidths = {
    0: 79,
    1:  76,
    2:  75}
    
    errorMessageTemplate =  " > > > {:<68}{:>18} < < < "
    leveledOutput = {0: " {:<80}{:>18} ",
                     1: "   {:<78}{:>18} ",
                     2: "     {:<76}{:>18} ",}
    
    def __init__(self, verbose=True, outputFile=None, suppressFromLevel=3):
        self.suppressLvl = suppressFromLevel
        self.verbose = verbose
    
    def message(self, message, senderIdentification, level=1):
        if level < self.suppressLvl:
            if self.verbose:
                print(self.leveledOutput[level].format(message, senderIdentification))
    
    def errorMessage(self, errorMessage, senderIdentification):
        print(self.errorMessageTemplate.format(errorMessage, senderIdentification))
        
    def printSeperationLine(self, ):
        if self.verbose:
            print('+'+'-'*98+'+')
            
    def printTable(self, table, senderIdentification, level=1, printHeaderRow=True):
        
        
        nCols = len(table[0])
        
        cellWidth = int( math.floor( self.outputWidths[level] / nCols - 1 ) )
        rowBar = '+' +  ( ('-' * cellWidth) + '+' ) * nCols
        rowString = ('|{:}'.format( '{:'+str(cellWidth)+'}' ) ) * nCols + '|'
        
        if printHeaderRow:
            self.message(rowBar, senderIdentification, level)
        for row in table:
            self.message(rowString.format(*row), senderIdentification, level)
        
        self.message(rowBar, senderIdentification, level)
        
    
    def setVerbose(self,):
        self.suppressLvl = 3
        
    def squelch(self, level):
        self.suppressLvl = level
    