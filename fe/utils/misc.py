    # -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 17:36:07 2016

@author: matthias
"""

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
        val = parts[1] if len(parts)>1 else 'True'
        resultDict[opt] = val
    return resultDict
    
    