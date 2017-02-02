import numpy as np
from fe.materials.umatlibrary cimport umatType, getUmat
cimport numpy as np
from libcpp.string cimport string


cdef public bint notificationToMSG(const string cppString):
    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
    print(cppString.decode('UTF-8'))
    return False

cdef extern from "uelCPS8RNonLocalSimpleUmat.h":
    void uelCPS8RNonLocalSimpleUmat(
            double rightHandSide[],           
            double KMatrix[],                            
            double stateVars[],                                         
            const int &nStateVars, 
            const double properties[],
            const int &nProperties,
            const double coordinates[],                 
            const double U_[],                                  
            const double dU_[],       
            const double time[2],                                       
            const double &dTime,                                        
            const int& elementNumber,                  
            double &pNewdT,         
            const int integerProperties[],
            const int &nIntegerProperties, 
            umatType umat,
            int nStateVarsUmat)
    
cdef void callUel(
        double[::1,:] Ke,
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
        umatType umat,
        int nStateVarsUmat):
    
    uelCPS8RNonLocalSimpleUmat(    &Pe[0],           
            &Ke[0][0],                            
            &stateVars[0],                                         
            stateVars.shape[0], 
            &properties[0],
            properties.shape[0],
            &coordinates[0],                 
            &UNew[0],                                  
            &dU[0],       
            &time[0],                                       
            dTime,                                        
            elNumber,                  
            pNewdT[0],         
            &intProperties[0],
            intProperties.shape[0], 
            umat,
            nStateVarsUmat
            )

cdef class Element:
    fields =                [["mechanical", "nonlocal damage"],
                             ["mechanical", "nonlocal damage"],
                             ["mechanical", "nonlocal damage"],
                             ["mechanical", "nonlocal damage"],
                             ["mechanical", "nonlocal damage"],
                             ["mechanical", "nonlocal damage"],
                             ["mechanical", "nonlocal damage"],
                             ["mechanical", "nonlocal damage"],                              
                             ] # fields identical for each node
    
    nNodes =                8
    nGaussPt =              4
    nDofPerEl =             24
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation  = np.array([0,1,3,4,6,7,9,10,12,13,15,16,18,19,21,22] + [2,5,8,11,14,17,20,23], dtype=int)  
    ensightType =           "quad8"
    
    cdef public nodes, nodeCoordinates, intProperties, elNumber
    
    cdef np.ndarray uelProperties, 
    cdef np.ndarray stateVars, stateVarsTemp
    cdef umatType umat
    cdef int nStateVars, nStateVarsUmat
    
    def __init__(self, nodes, elNumber):
        self.nodes = nodes
        self.nodeCoordinates = np.concatenate([ node.coordinates for node in nodes])
        self.intProperties = np.zeros(1, dtype=np.intc)
        self.elNumber = elNumber
    
    def setProperties(self, uelProperties, umatName, nStateVarsUmat):
        self.uelProperties = uelProperties
        self.nStateVarsUmat = nStateVarsUmat
        self.nStateVars = self.nGaussPt * (nStateVarsUmat + 13)
        self.stateVars = np.zeros(self.nStateVars)
        self.umat = getUmat(umatName.lower())
        
    def computeYourself(self, Ke, Pe, U, dU, time, dTime, pNewdT):
        self.stateVarsTemp = np.copy(self.stateVars)
        callUel(Ke,Pe,
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
    
    def acceptLastState(self,):
        self.stateVars = self.stateVarsTemp
    def resetToLastValidState(self,):
        pass
    
    resultIndices = {'stress': lambda nStateVarsUmat,kw, : np.arange(6) + (kw['location'])*nStateVarsUmat ,
                    'strain': lambda nStateVarsUmat,kw, : np.arange(6) + (kw['location'])*nStateVarsUmat +6,
                    'sdv':    lambda nStateVarsUmat,kw, : kw['indices'] +((kw['location'])-1)*nStateVarsUmat }  
                    
    def getResult(self, **kw):    
        stateVarIndices = self.resultIndices[kw['result']](self.nStateVarsUmat, kw)
        return self.stateVars[stateVarIndices]
