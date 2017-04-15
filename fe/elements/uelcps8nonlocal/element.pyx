import numpy as np
from fe.materials.umatlibrary cimport umatType, getUmat
cimport numpy as np
from libcpp.string cimport string

cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return False

cdef extern from "uelCPS8NonLocalSimpleUmat.h" nogil:
    void uelCPS8NonLocalSimpleUmat(
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
            const int &elementNumber,                  
            double& pNewdT,         
            const int integerProperties[],
            const int &nIntegerProperties, 
            umatType umat,
            int nStateVarsUmat)
    
cdef void callUel(
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
        double& dTime,
        int& elNumber,
        umatType umat,
        int nStateVarsUmat) nogil:
    
    uelCPS8NonLocalSimpleUmat(    &Pe[0],           
            &Ke[0],                            
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
                             ["mechanical", "nonlocal damage"],] # fields identical for each node
    
    nNodes =                8
    nGaussPt =              9
    nDofPerEl =             24
    sizeKe =                nDofPerEl * nDofPerEl
    dofIndicesPermutation = np.array([0,1,3,4,6,7,9,10,12,13,15,16,18,19,21,22] + [2,5,8,11,14,17,20,23], dtype=int)  
    ensightType =           "quad8"
    
    cdef public nodes, 
    cdef public int elNumber
    
    cdef double[::1] uelProperties, stateVars, stateVarsTemp, nodeCoordinates
    cdef umatType umat
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
        self.stateVars[:] = self.stateVarsTemp
        
    def resetToLastValidState(self,):
        pass
    
    resultIndices = {
                    'sdv':    lambda nStateVarsUmat,kw : kw['indices'] +(int(kw['location'])-1)*(nStateVarsUmat+12),
                    'stress': lambda nStateVarsUmat,kw : np.arange(6) + (int(kw['location'])-1)*(nStateVarsUmat+12) + nStateVarsUmat,
                    'strain': lambda nStateVarsUmat,kw : np.arange(6) + (int(kw['location'])-1)*(nStateVarsUmat+12) + nStateVarsUmat +6,
                    'all':     lambda nStateVarsUmat,kw: np.arange(4*(nStateVarsUmat+12))}  
                    
    def getResult(self, **kw):    
        stateVarIndices = self.resultIndices[kw['result']](self.nStateVarsUmat, kw)
        return np.asarray(self.stateVars)[stateVarIndices]
