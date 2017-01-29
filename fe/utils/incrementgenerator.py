# -*- coding: utf-8 -*-
"""
Created on Sat Jan  21 12:18:10 2017

@author: Matthias
"""

class IncrementGenerator:
    """
    Version 2 of the increment generator, 
    based on Version 1 of the  Abaqus-Workbench.
    Implementation as generator class !"""
    
    identification = "IncGen"
    
    def __init__(self, currentTime, stepLength, maxIncrement, minIncrement, maxNumberIncrements, journal):
        
        self.nPassedGoodIncrements = int(0)
        self.totalIncrements = int(0)
        self.maxIncrement = maxIncrement
        self.minIncrement = minIncrement
        self.maxNumberIncrements = maxNumberIncrements
        
        self.finishedStepProgress = 0.0
        self.increment = maxIncrement
        
        self.currentTime =      currentTime
        self.stepLength =       stepLength
        self.dT = 0.0
        self.journal = journal
        
    def generateIncrement(self):
        """ Returns the hexatuple consisting of 
            (increment number, increment fraction, finished step progress, 
            dT, increment start time of step, increment start time total)"""
        
        #zero increment
        yield (0, 0.0, 0.0, 0.0, 0.0, self.currentTime )
        
        while self.finishedStepProgress < (1.0-1e-15) and self.totalIncrements < self.maxNumberIncrements:
        
            if(self.nPassedGoodIncrements >= 3):
                self.increment *= 1.5
                if self.increment > self.maxIncrement:
                    self.increment = self.maxIncrement
                
            remainder = 1.0 - self.finishedStepProgress 
            if remainder < self.increment:
                self.increment = remainder
                
            dT =                            self.stepLength * self.increment
            startTimeOfIncrementInStep =    self.stepLength * self.finishedStepProgress     
            startTimeOfIncrementInTotal =   self.currentTime + startTimeOfIncrementInStep           
            self.finishedStepProgress +=    self.increment        
            
            self.totalIncrements +=         1
            self.nPassedGoodIncrements +=   1     
            
            yield (self.totalIncrements, 
                   self.increment, 
                   self.finishedStepProgress, 
                   dT, 
                   startTimeOfIncrementInStep, 
                   startTimeOfIncrementInTotal)
        
    def discardAndChangeIncrement(self, scaleFactor):
        """ Change increment size between minIncrement and 
        maxIncrement by a given scale factor."""
        
        if self.increment == self.minIncrement:
            self.journal.errorMessage("Cannot reduce increment size", self.identification)
            self.totalIncrements = self.maxNumberIncrements + 1
            return 
            
        self.finishedStepProgress -=    self.increment
        newIncrement = self.increment * scaleFactor
        if newIncrement > self.maxIncrement:
            self.increment = self.maxIncrement
        elif newIncrement < self.minIncrement:
            self.increment = self.minIncrement
        else:
            self.increment = newIncrement
            
        self.journal.message("Cutback to increment size {:}".format(self.increment),
                             self.identification, 2)
        self.nPassedGoodIncrements = 0
