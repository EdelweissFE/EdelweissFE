#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 20:26:22 2017

@author: matthias
"""
from fe.config.phenomena import fieldCorrectionTolerance, fluxResidualTolerance, fluxResidualToleranceAlternative
from fe.utils.misc import stringDict

def loadConfiguration(jobInfo):
    jobInfo['fieldCorrectionTolerance'] =           fieldCorrectionTolerance
    jobInfo['fluxResidualTolerance'] =              fluxResidualTolerance
    jobInfo['fluxResidualToleranceAlternative'] =   fluxResidualToleranceAlternative
    return jobInfo

def updateConfiguration(newConfiguration, jobInfo):
    target = newConfiguration['configuration']
    settings = stringDict( [ setting for line in newConfiguration['data'] for setting in line]  )
    jobInfo[target].update( { key : type(jobInfo[target][key])( val) for key, val in settings.items()  } )
