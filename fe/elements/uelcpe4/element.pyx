import numpy as np
from fe.materials.umatlibrary cimport pUmatType, getUmat
from fe.config.ueltypedefs cimport pSimpleUelWithUmatType
cimport numpy as np
from libcpp.string cimport string


cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return False
              
cdef extern from "userLibrary.h" namespace "userLibrary" nogil:
    pSimpleUelWithUmatType getSimpleUelWithUmatById(int id)

cdef class Element:
    fields =                [["mechanical"],
                             ["mechanical"],
                             ["mechanical"],
                             ["mechanical"],] # fields identical for each node
    nNodes =                4
    nGaussPt =              4
    nDofPerEl =             8
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation  = np.arange(0, 8, 1)
    ensightType =           "quad4"
    uelIdentification =     407
    
    cdef public nodes, 
    cdef public int elNumber
    
    cdef double[::1] uelProperties, stateVars, stateVarsTemp, nodeCoordinates
    cdef pUmatType umat
    cdef pSimpleUelWithUmatType uel
    cdef int nStateVars, nStateVarsUmat
    
    cdef int[::1] intProperties
    
    def __init__(self, nodes, elNumber):
        self.nodes = nodes
        self.nodeCoordinates = np.concatenate([ node.coordinates for node in nodes])
        self.elNumber = elNumber
        
    def setProperties(self, uelProperties, umatName, nStateVarsUmat):
        self.uelProperties = uelProperties
        self.nStateVarsUmat = nStateVarsUmat
        self.nStateVars = self.nGaussPt * (nStateVarsUmat + 12)
        self.stateVars = np.zeros(self.nStateVars)
        self.stateVarsTemp = np.zeros(self.nStateVars)
        self.umat = getUmat(umatName.lower())
        self.uel = getSimpleUelWithUmatById(self.uelIdentification)
        self.intProperties = np.empty(0, dtype=np.intc)#self.intProperties
        
    def computeYourself(self, 
                         double[::1] Ke, 
                         double[::1] Pe, 
                         const double[::1] U, 
                         const double[::1] dU, 
                         const double[::1] time, 
                         double dTime, 
                         double[::1] pNewdT):
        
        with nogil: # release the gil for parallel computing
            self.stateVarsTemp[:] = self.stateVars
            self.callUel(Ke,Pe,
                self.stateVarsTemp,
                U,
                dU,
                self.nodeCoordinates,
                self.uelProperties,
                self.intProperties,
                pNewdT,
                time,
                dTime,
                self.elNumber,
                self.umat,
                self.nStateVarsUmat)

    cdef void callUel(self,     
        double[::1] Ke,
        double[::1] Pe,
        double[::1] stateVars,
        double[::1] UNew,
        double[::1] dU,
        double[::1] coordinates,
        double[::1] properties,
        int[::1] intProperties,
        double[::1] pNewdT,
        double[::1] time,
        double dTime,
        int elNumber,
        pUmatType umat,
        int nStateVarsUmat) nogil:
        
            self.uel(&Pe[0], &Ke[0], &stateVars[0], stateVars.shape[0], &properties[0], 
                     properties.shape[0], &coordinates[0], &UNew[0], &dU[0], &time[0], 
                     dTime, elNumber, pNewdT[0], &intProperties[0], intProperties.shape[0], 
                     umat, nStateVarsUmat)
            
    def acceptLastState(self,):
        self.stateVars[:] = self.stateVarsTemp
        
    def resetToLastValidState(self,):
        pass
    
    resultIndices = {'stress': lambda nStateVarsUmat,kw, : np.arange(6) + (kw['location'])*nStateVarsUmat ,
                    'strain': lambda nStateVarsUmat,kw, : np.arange(6) + (kw['location'])*nStateVarsUmat +6,
                    'sdv':    lambda nStateVarsUmat,kw, : kw['indices'] +((kw['location'])-1)*nStateVarsUmat }  
                    
    def getResult(self, **kw):    
        stateVarIndices = self.resultIndices[kw['result']](self.nStateVarsUmat, kw)
        return np.asarray(self.stateVars)[stateVarIndices]
