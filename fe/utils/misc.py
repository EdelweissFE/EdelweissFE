    # -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 17:36:07 2016

@author: matthias
"""

import numpy as np

def flagDict(configLine):
    parts = [x.strip() for x in configLine.split("=")]
    opt = parts[0]
    val = True if (len(parts)>1 and parts[1] == "True") else False
    return {opt: val}
    
def stringDict(listOfStringAssigments):
    resultDict = {}
    
    for entry in listOfStringAssigments:
        parts = [x.strip() for x in entry.split("=")]
        opt = parts[0]
        val = '='.join(parts[1:]) if len(parts)>1 else 'True'
        resultDict[opt] = val
    return resultDict

def strToSlice(string):
    
    if ':' in string:
        idcs=[int (i) for i in string.split(':')]
        a, b = idcs
        return slice(a,b)
    else:
        return slice ( int (string), int (string)+1)

def isInteger(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
    
def filterByJobName(canditates, jobName):
    return [cand for cand in canditates if cand.get('jobName', 'defaultJob') == jobName]    


def mergeNumpyDataLines(multiLineData):
    flattenedMatProps = [p for row in multiLineData for p in row]
    return np.array(flattenedMatProps, dtype=np.float)