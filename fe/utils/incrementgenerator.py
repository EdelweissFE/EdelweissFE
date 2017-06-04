# -*- coding: utf-8 -*-
"""
Created on Sat Jan  21 12:18:10 2017

@author: Matthias
"""
from fe.utils.exceptions import ReachedMaxIncrements, ReachedMinIncrementSize

class IncrementGenerator:
    """ Version 2 of the increment generator, 
    based on Version 1 of the  Abaqus-Workbench.
    Implementation as generator class !"""
    
    identification = "IncGen"
    
    def __init__(self, currentTime, stepLength, maxIncrement, minIncrement, maxNumberIncrements, journal):
        
        self.nPassedGoodIncrements =    int(0)
        self.totalIncrements =          int(0)
        self.maxIncrement =             maxIncrement
        self.minIncrement =             minIncrement
        self.maxNumberIncrements =      maxNumberIncrements
        
        self.finishedStepProgress =     0.0
        self.increment =                maxIncrement
        self.allowedToIncreasedNext =   True
        
        self.currentTime =              currentTime
        self.stepLength =               stepLength
        self.dT =                       0.0
        self.journal =                  journal
        
    def generateIncrement(self):
        """ Returns the hexatuple consisting of 
            (increment number, increment fraction, finished step progress, 
            dT, increment start time of step, increment start time total)"""
        
        #zero increment; return value for first function call
        yield (0, 0.0, 0.0, 0.0, 0.0, self.currentTime )
        
        while self.finishedStepProgress < (1.0-1e-15):
        
            if ( self.totalIncrements >= self.maxNumberIncrements ):
                self.journal.errorMessage("Reached maximum number of increments", self.identification)
                raise ReachedMaxIncrements()
            
            if(self.nPassedGoodIncrements >= 3) and self.allowedToIncreasedNext:
                self.increment *= 1.5
                if self.increment > self.maxIncrement:
                    self.increment = self.maxIncrement
            self.allowedToIncreasedNext = True
                
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
            
    def preventIncrementIncrease(self,):
        """ May be called before an increment is requested, to prevent from 
        automatically increasing, e.g. in case of bad convergency """
        self.allowedToIncreasedNext = False
        
    def discardAndChangeIncrement(self, scaleFactor):
        """ Change increment size between minIncrement and 
        maxIncrement by a given scale factor."""
        
        if self.increment == self.minIncrement:
            self.journal.errorMessage("Cannot reduce increment size", self.identification)
            raise ReachedMinIncrementSize()
            
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
        self.totalIncrements -= 1
        self.nPassedGoodIncrements = 0
