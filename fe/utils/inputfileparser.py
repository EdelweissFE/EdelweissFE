#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""
import numpy as np
from os.path import dirname, join
import textwrap

dTypes = {int : "integer",
          float: "float",
          str: "string",          
          "numpy": "numpy array",
          "numpy int": "numpy int array",
          }
    
typeMappings = {    '*element':         ("definition of element(s)",
                        {   'type':     (str, "assign one of the types definied in the elementlibrary"),
                            'data':     ("numpy int", "Abaqus like element definiton lines"), }),

                    '*elSet':           ("definition of an element set",
                        {   'elSet':    (str, "name"),
                            'generate': (str, "set True to generate from data line 1: start-element, end-element, step"),
                            'data':     (int, "Abaqus like element set definiton lines"), }),

                    '*node':            ("definition of nodes",
                        {   
                            'data':     ("numpy", "Abaqus like node definiton lines: label, x, [y], [z]"), }),

                    '*nSet':            ("definition of an element set",
                        {   'nSet':     (str, "name"),
                            'generate': (str, "set True to generate from data line 1: start-node, end-node, step"),
                            'data':     (int, "Abaqus like node set definiton lines"), }),

                    '*section':         ("definition of an section",
                        {   'name':     (str, "name"),
                            'thickness':(float, "associated element set"),
                            'material': (str, "associated id of defined material"), 
                            'data':     (str, "list of associated element sets")}),

                    '*material':        ("definition of a material",
                        { 'name':       (str, "name of the property"),
                         'id':          (str, "name of the property"),
                         'statevars':   (int, "number of statevars"),
                         'data':        ("numpy", "material properties, multiline possible")}),
                                         
                     '*output':        ("define an output module",
                        { 
                         'jobName':     (str, "(optional), name of job, standard=defaultJob"),
                         'type':        (str, "output module "),
                         'data':        (str, "defintions lines for the output module")}),
                                         
                    '*job':                 ("definition of an analysis job",
                        {'name':            (str, "(optional) name of job, standard = defaultJob"),
                         'solver':          (str, "(optional) define solver, standard = NIST"),
                         'domain':          (str, "define spatial domain: 1d, 2d, 3d"),
                         'startTime':       (float,"(optional) start time of job, standard = 0.0"), 
                         }),
                                         
                     '*step':               ("definition of *job steps", 
                        {'stepLength':      (float, "time period of step"),
                         'jobName':         (str, "(optional), name of job, standard=defaultJob"),
                         'maxInc':          (float, "maximum size of increment"),
                         'minInc':          (float, "minimum size of increment"),
                         'maxNumInc':       (int, "maximum number of increments"),
                         'maxIter':         (int, "maximum number of increments"),
                         'data':            (str, "define step actions, which are handled by the corresponding stepaction modules") }),
                                         
                    '*include': ("(optional) load extra .inp file (fragment), use relative path to current .inp",
                        {'input':           (str, "filename")}),
                        
                }
                
def getMapType(kw, varName):
    kwDict = typeMappings.get(kw, (None,{}) )[1]    
    mType = kwDict.get(varName, [str])[0]
    return mType
    

def parseInputFile(fileName, currentKeyword = None, existingFileDict = None):
    """ Parse an Abaqus like input file to generate an dictionary with its content """
    
    if not existingFileDict:
        fileDict = { key : [] for key in typeMappings.keys()}
    else:
        fileDict = existingFileDict
    keyword = currentKeyword
    with open(fileName) as f:
        for l in f:
            lineElements = [x.strip() for x in l.split(",")]
            lineElements=list(filter(None,lineElements))
            if not lineElements or lineElements[0].startswith("**"):
                # line is a comment
                pass
            elif lineElements[0].startswith("*"):
                # line is keywordline
                lastkeyword = keyword
                keyword = lineElements[0]
                optionAssignments = lineElements[1:]
                
                objectentry = {}
                objectentry['data'] = []
                objectentry['inputFile'] = fileName # save also the filename of the original inputfile!
                    
                for ass in optionAssignments:
                    opts = ass.split("=")
                    optKey = opts[0].strip()
                    val = opts[1].strip()
                    mType = getMapType(keyword, optKey)
                    if mType is not None:
                        objectentry[optKey] = mType(val)
                    else:
                        objectentry[optKey] = val
                        
                #special treatment for *include:
                if keyword == "*include":
                    includeFile = objectentry['input']
                    parseInputFile(join(dirname(fileName), 
                                                          includeFile), 
                                                     currentKeyword=lastkeyword,
                                                     existingFileDict=fileDict)
                    keyword = lastkeyword
                else:
                    fileDict[keyword].append(objectentry)
                
            else:
                # line is a dataline
                data = lineElements
                mType = getMapType(keyword, "data")
                if mType is not None:
                    if mType == "numpy":
                        data = np.array([x for x in data], dtype = np.double)
                    elif mType == "numpy int":
                        data = np.array([x for x in data], dtype = np.int)
                    else:    
                        data = [mType(d) for d in data]
                fileDict[keyword][-1]['data'].append(data)
    
    # 'post processing' of the assembled filedict
#    for key, (desc, optionsDict) in typeMappings.items():
        # merge datalines, which are flags to a dictionary
#        if optionsDict['data'][0] == flagDict if 'data' in optionsDict else False:
#            for entry in fileDict[key]:
#                newDic  = {}
#                for l in entry['data']:
#                    for d in l:
#                        newDic.update(d)
#                entry['data'] = [[newDic]]
                
    return fileDict

def printKeywords():
    """ print the input file language set"""
    kwString = "    {:}    "
    kwDataString = "        {:22}{:20}"    
    
    wrapper = textwrap.TextWrapper(width=80,replace_whitespace=False)
    for kw, (kwDoc,optiondict) in sorted(typeMappings.items()):
        wrapper.initial_indent = kwString.format(str(kw))
        wrapper.subsequent_indent = " "*len(wrapper.initial_indent)
        print(wrapper.fill(kwDoc))
        print('')
        
        for key in sorted(optiondict.keys()):
            optionName = key
            dType, description = optiondict[key]
            wrapper.initial_indent = kwDataString.format(str(optionName),dTypes[dType])
            wrapper.subsequent_indent = " "*len(wrapper.initial_indent)
            print(wrapper.fill(description))
        print("\n")
    
def mergeDictDataLines(multiLineData):
    d = {key:val for (key, val) in multiLineData}
    return d
    
