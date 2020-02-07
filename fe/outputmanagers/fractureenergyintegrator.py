#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 08:26:01 2019

@author: Matthias Neuner

A simple integrator to compute the fracture energy

Datalines:
"""
documentation = {
                'forceFieldOutput' : 'fieldOutput for force (with time history)',
                'displacementFieldOutput' : 'fieldOutput for displacement (with time history)',
                 'fractureArea' : '(math expression for) area of fracture',
                 }

from fe.outputmanagers.outputmanagerbase import OutputManagerBase
from fe.utils.misc import stringDict
from fe.utils.math import createMathExpression

import numpy as np

class OutputManager(OutputManagerBase):
    """ Simple Integrator for fracture energies"""
    
    identification = "FEI"
    printTemplate = "{:}, {:}: {:}"
    
    def __init__(self, name, definitionLines, jobInfo, modelInfo, fieldOutputController, journal, plotter):
        self.journal = journal
        self.monitorJobs = []
        self.fieldOutputController = fieldOutputController

        mergedDefinition = [x for l in definitionLines for x in l ]
        
        defDict = stringDict(mergedDefinition)
        self.fpF = fieldOutputController.fieldOutputs [ defDict['forceFieldOutput'] ]
        self.fpU = fieldOutputController.fieldOutputs [ defDict['displacementFieldOutput'] ]
        self.A = createMathExpression ( defDict['fractureArea'] )( 0.0 )
        self.fractureEnergy = 0.0

    def initializeStep(self, step, stepActions, stepOptions):
        pass
    
    def finalizeIncrement(self, U, P, increment):
        pass
    
    def finalizeStep(self, U, P):
        pass
    
    def finalizeJob(self,U, P):
        self.fractureEnergy = np.trapz( self.fpF.getResultHistory(), x = self.fpU.getResultHistory() ) / self.A
        self.journal.message("integrated fracture energy: {:3.4f}".format ( self.fractureEnergy), self.identification)
