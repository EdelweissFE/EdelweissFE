#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  ---------------------------------------------------------------------
#
#  _____    _      _              _         _____ _____ 
# | ____|__| | ___| |_      _____(_)___ ___|  ___| ____|
# |  _| / _` |/ _ \ \ \ /\ / / _ \ / __/ __| |_  |  _|  
# | |__| (_| |  __/ |\ V  V /  __/ \__ \__ \  _| | |___ 
# |_____\__,_|\___|_| \_/\_/ \___|_|___/___/_|   |_____|
#                                                       
# 
#  Unit of Strength of Materials and Structural Analysis
#  University of Innsbruck,
#  2017 - today
# 
#  Matthias Neuner matthias.neuner@uibk.ac.at
# 
#  This file is part of EdelweissFE.
# 
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
# 
#  The full text of the license can be found in the file LICENSE.md at
#  the top level directory of EdelweissFE.
#  ---------------------------------------------------------------------
"""
Created on Mon Apr 10 20:26:22 2017

@author: Matthias Neuner
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
