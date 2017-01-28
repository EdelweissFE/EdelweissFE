#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 11:30:59 2017

@author: matthias
"""

class Journal:
    """
    This class provides an interface to present messages to the user via console
    output and/or file output.
    Information messages can be sorted by the importance level. 
    Suppressing certain levels of output is possible.
    """
    digitsLevel0Output = 80
    digitsLevel1Output = 78
    digitsLevel2Output = 76
    errorMessageTemplate =  " > > > {:<68}{:>18} < < < "
    leveledOutput = {
                     0: " {:<80}{:>18} ",
                     1: "   {:<78}{:>18} ",
                     2:"     {:<76}{:>18} ",
    }
    def __init__(self, outputFile=None, suppressFromLevel=3):
        self.suppressLvl = suppressFromLevel
    
    def message(self, message, senderIdentification, level=1):
        if level < self.suppressLvl:
            print(self.leveledOutput[level].format(message, senderIdentification))
    
    def errorMessage(self, errorMessage, senderIdentification):
        print(self.errorMessageTemplate.format(errorMessage, senderIdentification))
        
    def printSeperationLine(self, ):
        print('+'+'-'*98+'+')
    
    def setVerbose(self,):
        self.suppressLvl = 3
        
    def squelch(self, level):
        self.suppressLvl = level
    