#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 14:22:48 2017

@author: matthias
"""
from fe.outputmanagers.outputmanagerbase import OutputManagerBase

import os
import datetime
import numpy as np
from collections import defaultdict, OrderedDict
from fe.utils.misc import stringDict
import fe.config.phenomena

def writeCFloat(f, ndarray):
    np.asarray(ndarray, dtype=np.float32).tofile(f)
def writeCInt(f, ndarray):
    np.asarray(ndarray, dtype=np.int32).tofile(f)
def writeC80(f, string):
    np.asarray(string, dtype='a80').tofile(f)

variableTypes = {"scalar" : 1,
                 "vector" : 3,
                 "tensor": 9}

ensightPerNodeVariableTypes = {
                               1 : 'scalar per node',
                               3 : 'vector per node',
                               6 : 'tensor per node',
                               9 : 'tensor9 per node'}  
                                   
ensightPerElementVariableTypes = {
                                  1 : 'scalar per element',
                                  3 : 'vector per element',
                                  6 : 'tensor per element',
                                  9 : 'tensor9 per element'}
                 
class EnsightUnstructuredPart:
    """ define an unstructured part, by a list of nodes and a dictionary of elements.
    Each dictionary entry consists of a list of tuples of elementlabel and nodelist: 
    {strElementType : [ ( intLabel = None, [nodeList] ) ]}"""
    def __init__(self, description, partNumber, nodes, elementTree, ):
        
        self.structureType = "coordinates"
        self.elementTree = elementTree 
        
        self.nodes = nodes
        self.nodeLabels =  np.asarray([node.label for node in nodes ], np.int32 )
            
            
        self.description = description  #string, describing the part; max. 80 characters
        self.partNumber = partNumber 
        self.nodeCoordinateArray = np.asarray( [node.coordinates for node in nodes] ) 
        
    def writeToFile(self, binaryFileHandle, printNodeLabels=True, printElementLabels=True):
        
        if len(self.nodeCoordinateArray.shape) > 1  and self.nodeCoordinateArray.shape[1] < 3:
            extendTo3D = True
        else:
            extendTo3D = False
            
        nNodes = self.nodeCoordinateArray.shape[0]
        f = binaryFileHandle  
        
        writeC80(f, 'part')
        writeCInt(f, self.partNumber)
        writeC80(f, self.description)
        writeC80(f, 'coordinates')
        writeCInt(f, nNodes)
        
        if printNodeLabels and self.nodeLabels is not None:                 
            writeCInt(f, self.nodeLabels)
            
        writeCFloat(f, self.nodeCoordinateArray.T)
        
        if extendTo3D:
            writeCFloat(f, np.zeros(self.nodeCoordinateArray.shape[0] * (3 - self.nodeCoordinateArray.shape[1])))
                
        for elemType, elements in self.elementTree.items():
            writeC80(f, elemType)
            writeCInt(f, len(elements))
            if printElementLabels:
                writeCInt(f, np.asarray([element.elNumber for element in elements.keys()], np.int32))
                
            for nodeIndices in elements.values():             
                writeCInt(f, np.asarray(nodeIndices, np.int32)[:]+1 )

#class EnsightPointPart(EnsightUnstructuredPart):
#    """derived UnstructuredPart to represent a single point"""
#    def __init__(self, partNumber, coordinates, description="single_point"):
#        super().__init__(description, partNumber, {"point" : [(1, [0])] }, [ {'label': 1, 'coords':coordinates, } ])

class EnsightTimeSet:
    """ defines a set which may be used by EnsightGeometry, EnsightStructuredPart, EnsightUnstructuredPart and is written into the case file"""
    def __init__(self,number=1, description="timeStepDesc",  fileNameStartNumber=0, fileNameNumberIncrement=1, timeValues = None):
        self.number = number
        self.description = description
        self.fileNameStartNumber = fileNameStartNumber
        self.fileNameNumberIncrement = fileNameNumberIncrement
        self.timeValues = timeValues if timeValues is not None else []
        
class EnsightGeometryTrend:
    """ container class for the time dependent evolution of the geometry,
        establishes the connection between the geometry entities and a EnsightTimeSet"""
    def __init__(self, ensightTimeSet, ensightGeometryList = None):
        self.timeSet = ensightTimeSet
        self.geometryList = ensightGeometryList if ensightGeometryList is not None else []
        
class EnsightGeometry:
    """ container class for one or more EnsightParts at a certain time state, handles also the file writing operation"""
    def __init__(self, name="geometry", descriptionLine1 ="", descriptionLine2 ="", ensightPartList = None , nodeIdOption="given", elementIdOption="given"):
        self.name = name
        self.descLine1 = descriptionLine1
        self.descLine2 = descriptionLine2
        self.partList = ensightPartList if ensightPartList is not None else []
        self.nodeIdOption = nodeIdOption
        self.elementIdOption = elementIdOption
    
    def writeToFile(self, fileHandle):
        f = fileHandle
        writeC80(f, self.descLine1)
        writeC80(f, self.descLine2)
        writeC80(f, "node id "+self.nodeIdOption)
        writeC80(f, "element id "+self.elementIdOption)
        
        if self.nodeIdOption == "given" or self.nodeIdOption == "ignore":
            printNodeLabels = True
        else:# assign or off
            printNodeLabels = False
            
        if self.elementIdOption == "given" or self.nodeIdOption == "ignore":
            printElementLabels = True
        else: # assign or off 
            printElementLabels = False
            
        for part in self.partList:
            part.writeToFile(f, printNodeLabels, printElementLabels)
        


class EnsightVariableTrend:
    """ container class for the time dependent evolution of one variable,
        establishes the connection between EnsightVariable entities and a EnsighTimeSet"""
    def __init__(self, ensightTimeSet, variableName,  ensightVariableList = None, variableType = "scalar per node", description = "variableTrendDescription"):
        self.timeSet = ensightTimeSet
        self.variableName = variableName
        self.variableList = ensightVariableList if ensightVariableList is not None else []
        self.variableType = variableType
        self.description = description
        
class EnsightPerNodeVariable:
    """ container class for data for one certain variable, defined for one or more parts (classification by partID), at a certain time state.
    For each part the structuretype ("coordinate" or "block") has to be defined.
    Each part-variable assignment is defined by a dictionary entry of type: { EnsightPart: np.array(variableValues) }"""
    def __init__(self, name, variableDimension, ensightPartsDict = None) :
        self.name = name
        self.description = name
        self.partsDict = ensightPartsDict or {} # { EnsightPart: np.array(variableValues) }
        self.variableDimension = variableDimension
        self.varType = ensightPerNodeVariableTypes[variableDimension]
        
    def writeToFile(self, fileHandle, ):
        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, (structureType,  values) in self.partsDict.items():          
            writeC80(f, 'part')
            writeCInt(f, ensightPartID)
            writeC80(f, structureType)
            writeCFloat(f, values.T)
            if values.shape[1] < self.variableDimension:
                writeCFloat(f, np.zeros( (values.shape[0],self.variableDimension - values.shape[1])))
            
class EnsightPerElementVariable:
    """ Container class for data for one certain variable, defined for one or more parts (classification by partID), at a certain time state.
    For each part the structuretype ("coordinate" or "block") has to be defined.
    Each part-variable assignment is defined by a dictionary entry of type: { EnsightPart: np.array(variableValues) }"""
    
    def __init__(self, name, variableDimension, ensightPartsDict = None, ) :
        self.name = name
        self.description = name
        self.partsDict = ensightPartsDict or {} 
        self.varType = ensightPerElementVariableTypes[variableDimension]
        
    def writeToFile(self, fileHandle):
        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, elTypeDict in self.partsDict.items():          
            writeC80(f, 'part')
            writeCInt(f, ensightPartID)
            for elType, values in elTypeDict.items():
                writeC80(f, elType)
                writeCFloat(f, values.T)

class EnsightChunkWiseCase:
    
    def __init__(self, caseName, directory = '', writeTransientSingleFiles = True):
        self.directory = directory
        self.caseName = caseName
        self.caseFileNamePrefix = os.path.join(directory, caseName)  
        self.writeTransientSingleFiles = writeTransientSingleFiles
        self.timeAndFileSets = {}
        self.geometryTrends = {}
        self.variableTrends = {}
        self.fileHandles = {}
        
        os.mkdir( self.caseFileNamePrefix )

    def setCurrentTime(self, timeAndFileSetNumber, timeValue):
        if not timeAndFileSetNumber in self.timeAndFileSets:
            self.timeAndFileSets[timeAndFileSetNumber] = EnsightTimeSet(timeAndFileSetNumber, 'noDesc', 0,1)
        tfSet = self.timeAndFileSets[timeAndFileSetNumber]
        tfSet.timeValues.append(timeValue)
        
    def writeGeometryTrendChunk(self, ensightGeometry, timeAndFileSetNumber=1):
        
        if ensightGeometry.name not in self.fileHandles:
            fileName = os.path.join (self.caseFileNamePrefix, ensightGeometry.name + ".geo",)
            
            self.fileHandles[ensightGeometry.name] = open(fileName, mode='wb')
        
        f = self.fileHandles[ensightGeometry.name]
        

        if not ensightGeometry.name in self.geometryTrends:
            self.geometryTrends[ensightGeometry.name] = timeAndFileSetNumber
            writeC80(f, 'C Binary')
                
        if self.writeTransientSingleFiles:
            writeC80(f, 'BEGIN TIME STEP')
            ensightGeometry.writeToFile(f)
            writeC80(f, 'END TIME STEP')
        
    def writeVariableTrendChunk(self, ensightVariable, timeAndFileSetNumber=2):
        
        if ensightVariable.name not in self.fileHandles:
            
            fileName = os.path.join(self.caseFileNamePrefix, ensightVariable.name + ".var")
            
            self.fileHandles[ensightVariable.name] = open(fileName, mode='wb')
        
        f = self.fileHandles[ensightVariable.name]

        
        if not ensightVariable.name in self.variableTrends:
            self.variableTrends[ensightVariable.name] = timeAndFileSetNumber, ensightVariable.varType
            writeC80(f, 'C Binary')
                
        if self.writeTransientSingleFiles:
            writeC80(f, 'BEGIN TIME STEP')
            ensightVariable.writeToFile(f)
            writeC80(f, 'END TIME STEP')
                    
    def finalize(self, replaceTimeValuesByEnumeration=True, closeFileHandes=True):
        caseFName = self.caseFileNamePrefix+'.case'
        
        if closeFileHandes:
            for f in self.fileHandles.values():
                f.close()
        
        with open(caseFName ,mode='w') as cf:
            cf.write("FORMAT\n")
            cf.write("type: ensight gold\n")
        
            cf.write("TIME\n")
            for setNum, timeSet in self.timeAndFileSets.items():
                cf.write("time set: "+str(setNum)+" noDesc\n")
                cf.write("number of steps: "+str(len (timeSet.timeValues))  + "\n")
                cf.write("filename start number: " + str(timeSet.fileNameStartNumber) +"\n")
                cf.write("filename increment: " +str(timeSet.fileNameNumberIncrement) +"\n")
                cf.write("time values: ")
                for i, timeVal in enumerate(timeSet.timeValues):
                    if not replaceTimeValuesByEnumeration:
                        cf.write('{:1.8e}'.format(timeVal) +"\n")
                    else:
                        cf.write('{:}'.format(i) +"\n")
                
            if self.writeTransientSingleFiles:
                cf.write("FILE\n")
                for timeSet in self.timeAndFileSets.values():
                    cf.write("file set: {:} \n".format(timeSet.number))
                    cf.write("number of steps: {:} \n".format(len(timeSet.timeValues)))
            
            cf.write("GEOMETRY\n")
            for geometryName, tAndFSetNum in self.geometryTrends.items():         
                cf.write("model: {:} \n".format(os.path.join ( self.caseFileNamePrefix, geometryName + ".geo" ) ))
                
            cf.write("VARIABLE\n")
            for variableName, (tAndFSetNum, variableType) in self.variableTrends.items():
                cf.write("{:}: {:} {:} {:} {:}.var\n".format(
                    variableType,
                    tAndFSetNum,
                    tAndFSetNum,
                    variableName,
                   os.path.join( self.caseFileNamePrefix, variableName)  ))
                
def createUnstructuredPartFromElementSet(setName, elementSet, partID):
    """ Determines the element and node list for an Ensightpart from an 
    element set. The reduced, unique node set is generated, as well as 
    the element to node index mapping for the ensight part."""
    
    nodeCounter = 0
    partNodes = OrderedDict() # node -> index in nodelist
    elementDict = defaultdict(OrderedDict)
    for element in elementSet:
        elNodeIndices = []
        for node in element.nodes:
            # if the node is already in the dict, get its index, 
            # else insert it, and get the current idx = counter. increase the counter
            idx = partNodes.setdefault(node, nodeCounter)
            elNodeIndices.append(idx)
            if idx == nodeCounter:
                # the node was just inserted, so increase the counter of inserted nodes
                nodeCounter+=1
        elementDict[element.ensightType][element] = elNodeIndices 
    
    return EnsightUnstructuredPart(setName, partID, partNodes.keys(), elementDict)

class OutputManager(OutputManagerBase):
    identification = "Ensight Export"
    
    def __init__(self, name, definitionLines, jobInfo, modelInfo, journal):
        self.name = name
        
        self.finishedSteps = 0
        self.intermediateSaveInterval = 0
        self.intermediateSaveIntervalCounter = 0
        self.domainSize = jobInfo['domainSize']
        self.journal = journal
        
        exportName = '{:}_{:}'.format(name,  datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))
        
        self.elSetToEnsightPartMappings = {}
        self.ensightCase = EnsightChunkWiseCase(exportName)
        
        self.perNodeJobs = []
        self.perElementJobs = []
                    
        elementSets = modelInfo['elementSets']
        
        elSetParts = []
        partCounter = 1
        for setName, elSet in elementSets.items():
            elSetPart = createUnstructuredPartFromElementSet(setName, elSet, partCounter)
            self.elSetToEnsightPartMappings[setName] = elSetPart
            elSetParts.append(elSetPart)
            partCounter += 1
            
        geometry = EnsightGeometry("geometry", "Edelweiss_FE", "*export*", ensightPartList=elSetParts)
    
        geometryTimesetNumber = None 
        self.ensightCase.writeGeometryTrendChunk(geometry, geometryTimesetNumber)
        
        for defLine in definitionLines:
            definition = stringDict(defLine)
            if 'create' in definition:
                varType = definition['create']
                
                if varType == 'perNode':
                    if 'elSet' in definition:
                        setName = definition['elSet'].strip()
                        perNodeJob = {}
                        perNodeJob['name'] = definition['name']
                        perNodeJob['part'] =  self.elSetToEnsightPartMappings[setName]
                        field = definition['field']
                        perNodeJob['result'] = definition['type']
                        perNodeJob['varSize'] = variableTypes[ fe.config.phenomena.phenomena[field] ]
                        perNodeJob['indices'] = np.asarray([node.fields[field] for node in perNodeJob['part'].nodes]).ravel()
                        perNodeJob['dimensions'] = fe.config.phenomena.getFieldSize(field, self.domainSize)
                        self.perNodeJobs.append(perNodeJob)
                        
                if varType == 'perElement':
                    perElementJob = {}
                    perElementJob['part'] =  self.elSetToEnsightPartMappings[definition['elSet']]
                    perElementJob['result'] = definition['result']
                    
                    if 'index' in definition:
                        idcs = definition['index']
                        if ':' in definition['index']:
                            idcs=[int (i) for i in idcs.split(':')]
                            perElementJob['idxStart'], perElementJob['idxStop'] = idcs
                        else:
                            idx = int (idcs )
                            perElementJob['idxStart'], perElementJob['idxStop'] = idx, idx+1
                        
                    if 'gaussPt' in definition:
                        perElementJob['gaussPt'] = int( definition['gaussPt'] )
                        
                    perElementJob['name'] = definition.get('name', perElementJob['result'])
                    
                    part = perElementJob['part']
                    perElementJob['permanentElResultMemory'] = {}
                    for ensElType, elements in part.elementTree.items():
                        perElementJob['permanentElResultMemory'][ensElType] = [el.getPermanentResultPtr(**perElementJob) for el in  elements.keys()]
                    
                    self.perElementJobs.append(perElementJob)
                    
    def initializeStep(self, step, stepActions):
        if 'EnsightOptions' in stepActions and self.name in stepActions['EnsightOptions']:
            options = stepActions['EnsightOptions'][self.name]
            self.intermediateSaveInterval = int(options.get('intermediateSaveInterval', 
                                                            self.intermediateSaveInterval))
            
    def finalizeIncrement(self, U, P, increment):
        incNumber, incrementSize, stepProgress, dT, stepTime, totalTime = increment
        
        if dT <= 0.0 and self.finishedSteps > 0:
            self.journal.message("skipping output for zero-increment in step {:}"
                                 .format(self.finishedSteps+1), self.identification, 1)
            return
        
        self.ensightCase.setCurrentTime(1, totalTime + dT)
        
        for perNodeJob in self.perNodeJobs:
            resultIndices = perNodeJob['indices']
            resultTypeLength = perNodeJob['varSize'] 
            resultDimensions = perNodeJob['dimensions']
            jobName = perNodeJob['name']
            nodalVarTable = U[resultIndices].reshape( (-1, resultDimensions)  )
            partsDict = {perNodeJob['part'].partNumber : ('coordinates', nodalVarTable)}
            enSightVar = EnsightPerNodeVariable(jobName, resultTypeLength, partsDict)
            self.ensightCase.writeVariableTrendChunk(enSightVar, 1)
            del enSightVar
            
        for perElementJob in self.perElementJobs:
            name = perElementJob['name']            
            part = perElementJob['part']
            varDict = {}
            for ensElType, elements in part.elementTree.items():
                varDict[ensElType] = np.asarray(perElementJob['permanentElResultMemory'][ensElType])
                if len(varDict[ensElType].shape)==1:
                    dimension=1
                else:
                    dimension = varDict[ensElType].shape[1]

            partsDict = {part.partNumber : varDict}
            enSightVar = EnsightPerElementVariable(name, dimension, partsDict)
            self.ensightCase.writeVariableTrendChunk(enSightVar, 1)
            del enSightVar  
        
        # intermediate save of the case
        if self.intermediateSaveInterval:
            if self.intermediateSaveIntervalCounter == self.intermediateSaveInterval:
                self.ensightCase.finalize(replaceTimeValuesByEnumeration=False, closeFileHandes=False)
                self.intermediateSaveIntervalCounter = 0
            self.intermediateSaveIntervalCounter +=1
            
    def finalizeStep(self, U, P,):
        self.finishedSteps += 1
        
    def finalizeJob(self, U, P,):
        self.ensightCase.finalize(replaceTimeValuesByEnumeration=False)
        
def printDocumentation():
    print("""Export a job to an Ensight Gold Case
data line: create=perNode|perElement
 - perNode: elSet: set, for which the per node job is created
            name: variable export name
            type: U|P (flow|effort)
            field: result field
 - perElement: elSet: set, for which the per element job is created
               name: variable export name
               result: element dependent result name (e.g. sdv, stress, strain),
                       passed to element
               index: index (or slice) for the result extraction, 
                      passed to element
               (gaussPt): optional, id of gaussPt, passed to element
               
EnsightOptions in step actions:
 - intermediateSaveInterval=N : intermediate save every N increments
 """)
