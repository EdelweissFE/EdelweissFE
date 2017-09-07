#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 16:29:29 2017

@author: c8441146
"""
import numpy as np
from scipy.optimize import curve_fit

from fe.fecore import finitElementSimulation
from fe.utils.misc import stringDict, mergeNumpyDataLines

import matplotlib.pyplot as plt

materialParametersLists = {
    'untereggerrockmass': ['E', 'nu', 'fcu', 'm0', 'e', 'mg1', 'fcy', 'Ah', 'Bh', 'Ch', 'Dh', 'Gh', 'As', 'Bs', 'GSI', 'D', 'omegaMax', 'NLradius', 'NLweight', 'epsilonF'],                       
#    'shotleon':[ 'q1', 'q2', 'q3', 'q4', 'nu', 'E1', 'nKelvin', 'kelvinMin', 'strainShrInf', 'TauShr', 'kHum', 'fcu28', 'fcu1', 'ratioFcy', 'ratioFbu', 'ratioFtu', 'eCPP1', 'eCPP8', 'eCPP24', 'Df', 'Ah', 'Ch', 'Dh', 'As', 'GFI', 'enableDamage', 'dTStatic', 'castTime', 'timeToDays',],
#    'meschke':['E28', 'nu', 'E1', 'fc28', 'fc1', 'fb_fc', 'fcy_fc', 'varEps', 'tE', 'deltaTE', 'Gfpython buildWorkbench.py build_ext -iI', 'viscosity', 'shrStrainInf', 'TauShr', 'kHum', 'timeToHours',],
#    'schaedlichschweiger':['E28', 'nu', 'fc28', 'E1', 'fc1', 'eCpp1', 'eCpp8', 'eCpp24', 'fc0n', 'fcfn', 'fcun', 'ftun', 'Gc28', 'Gt28', 't50Cr', 'strainInfS', 't50Shr', 'phiCr', 'phiMCMax', 'psiDilatancy', 'aDuct', 'varEps', 'varEpsR', 'eR', 'timeToHours'],
    }


def createAbaqusMaterialInputDataLines(materialName, matParameters):
    """ create commented, Abaqus-valid material input dataLines""" 
    dataLines = []
    parameterNames = materialParametersLists.get(materialName.lower(), False)
    if not parameterNames:
        parameterNames = ['Parameter #{:}'.format(i) for i in range(len(matParameters)) ]
    n = 8
    for i in range(0, len(matParameters), n):
        parNamesChunk = parameterNames[i:i+n]
        parChunk =      matParameters[i:i+n]
        formatString = '{:>15},'*len(parChunk)
        dataLines.append('**'+formatString.format(*parNamesChunk))
        dataLines.append(('  '+'{:>15.4e},'*len(parChunk)).format(*parChunk))
        
    return dataLines
    
def evaluateJob(inputFile, xFieldOutput, yFieldOutput, matParams, paramIndices, xData, *params): 
    """ wrapper to generate vectorial data from FEM simulation for given xData 
    and for certain material parameters to be used within the parameter fitting
    """
    for idx, param in zip(paramIndices, params):
        matParams[idx] = param
        
    success, U, P, fieldOutputController  = finitElementSimulation(inputFile, verbose=False)
           
    xValues = fieldOutputController.fieldOutputs[ xFieldOutput ].result 
    yValues = fieldOutputController.fieldOutputs[ yFieldOutput ].result 

    if success:
#    
#         flip if x is not ascending (necessary for interpolation by np.interp)    
        if xValues[0]>xValues[1]:
            xValues = xValues[::-1]
            yValues = yValues[::-1]
#        
        yValues = np.interp(xData, xValues, yValues)
#        
    else:   #FEM failed, return only zeros .. TODO: better error handling?
        yValues = np.zeros_like(xData)
#    
    return yValues

def parameterIdentification(inputFile):
    """ parameter identification for given x-y data based on the numpy curve_fit 
        function; 
        
        call parameteridentification by: 
                edelweissFE yourInpFile.inp --parameterIdentification True
                
    """
    for identificationJob in inputFile['*parameterIdentification']:
        
        for idx in range(len(inputFile['*material'])):
            if inputFile['*material'][idx]['id']==identificationJob['id']:
                material = inputFile['*material'][idx] 
            
        material['data'] = mergeNumpyDataLines(material['data'])
        
        xVals = np.loadtxt(identificationJob['xData'])
        yVals = np.loadtxt(identificationJob['yData'])     
        
        indices = []
        estimations = []
        lowerBounds = []
        upperBounds = []

        for dataline in identificationJob['data']:
            if 'type=x' in dataline:
                xFieldOutput = stringDict(dataline)['fieldOutput']
            elif 'type=y' in dataline:
                yFieldOutput = stringDict(dataline)['fieldOutput']
            else: 
                dataline = stringDict(dataline)            
                indices.append( int(dataline.get('idx')))
                estimations.append(float(dataline.get('start')))
                lowerBounds.append(float(dataline.get('min')))
                upperBounds.append(float(dataline.get('max')))

        matParams = material['data']
        
        inputFile['*output'] = []
        inputFile['*configurePlots'] = []

        callBackFunc = lambda xData, *params: evaluateJob(inputFile, xFieldOutput, yFieldOutput, matParams, indices, xData, *params)
        
        # sorting !
        sortedIndices = xVals.argsort()
        xVals = xVals[sortedIndices]
        yVals = yVals[sortedIndices]

        popt, pcov = curve_fit(callBackFunc, xVals, yVals, p0 = estimations, 
                       bounds=(lowerBounds, upperBounds), 
                        verbose=2,
                        gtol=1e-11,
                        xtol=1e-11,)

        for idx, param in zip(indices, popt):
            matParams[idx] = param      
            
        print('*'*80)
        print("new parameters:")      
        headerComment="**{:>8}{:>17}{:>17}{:>17}{:>17}"
        print(headerComment.format("#","initial value","lowerBound", "upper bound", "optimal value"))

        materialParameterNames = materialParametersLists.get(material['name'], False)
        if materialParameterNames:
            names = [materialParameterNames[idx] for idx in indices]
        else:
            names = indices

        parameterLine = ""    
        for idx, est, lower, upper, ident in zip(names, 
                                                         estimations, 
                                                         lowerBounds, 
                                                         upperBounds, 
                                                         popt):
            parameterLine += "**{:>8}{:>17}{:>17}{:>17}{:>17.4e}".format(idx, est, lower, upper, ident)
        
        print(parameterLine)
        np.set_printoptions(suppress=True, precision=5, )
        print('*'*80)
            
        if identificationJob['file']:
            with open(identificationJob['file'], 'w+' ) as f:
                f.write(parameterLine+'\n')
                parameterDataLines = createAbaqusMaterialInputDataLines(material['name'], matParams)
                for line in parameterDataLines:
                    f.write(line + '\n')
            
        if identificationJob.get('plot',''):
            plt.plot(xVals, yVals, 'r--', xVals, callBackFunc(xVals, np.asarray(estimations)), 'g:', xVals, callBackFunc(xVals, popt), 'b-.')
            plt.legend(['given data', 'initial val', 'fitted curve'])
            plt.grid()
            plt.show()

    return
