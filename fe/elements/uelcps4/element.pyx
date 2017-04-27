import numpy as np
from fe.materials.umatlibrary cimport pUmatType, getUmat
from fe.config.ueltypedefs cimport pSimpleUelWithUmatType
cimport numpy as np
from libcpp.string cimport string
from fe.elements.AbstractBaseElements.CPPBackendedElement cimport BackendedElement

cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return False

cdef extern from "NISTParallelizableBackendElement.h":
    cdef cppclass NISTParallelizableBackendElement:
        NISTParallelizableBackendElement(int elNumber, 
                                        const double* coordinates, 
                                        double *stateVars,
                                        int nStateVars, 
                                        const double* properties,
                                        int nProperties,
                                        const int* intProperties,
                                        int nIntProperties,
                                        pUmatType umat,
                                        int nStateVarsUmat,
                                        pSimpleUelWithUmatType uel)

        void computeYourself(double* Pe, double* Ke,  const double* UNew, const double* dU,  const double time[], double dTime, double &pNewDT )
        void acceptLastState()
              
cdef extern from "userLibrary.h" namespace "userLibrary" nogil:
    pSimpleUelWithUmatType getSimpleUelWithUmatById(int id)

cdef class Element(BackendedElement):
    fields =                [["mechanical"], # node 1
                             ["mechanical"], # node 2
                             ["mechanical"], # node 3
                             ["mechanical"],]# node 4
    nNodes =                4
    nGaussPt =              4
    nDofPerEl =             8
    sizeKe =                8 * 8
    dofIndicesPermutation = np.array([0,1,2,3,4,5,6,7])
    ensightType =           "quad4"
    uelIdentification =     402
    
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
        self.umat = getUmat(umatName.lower())
        self.uel = getSimpleUelWithUmatById(self.uelIdentification)
        self.intProperties = np.empty(0, dtype=np.intc)
        self.backendElement = new NISTParallelizableBackendElement(self.elNumber ,
                                                              &self.nodeCoordinates[0], 
                                                              &self.stateVars[0],
                                                              self.nStateVars,
                                                              &self.uelProperties[0], 
                                                              self.uelProperties.shape[0], 
                                                              &self.intProperties[0], 
                                                              self.intProperties.shape[0], 
                                                              self.umat, 
                                                              self.nStateVarsUmat, 
                                                              self.uel)
        
    def computeYourself(self, 
                         double[::1] Ke, 
                         double[::1] Pe, 
                         const double[::1] U, 
                         const double[::1] dU, 
                         const double[::1] time, 
                         double dTime, 
                         double[::1] pNewdT):
        self.backendElement.computeYourself(&Pe[0], &Ke[0], &U[0], &dU[0], &time[0], 
             dTime,  pNewdT[0])
    
    def acceptLastState(self,):
        self.backendElement.acceptLastState()
        
    def resetToLastValidState(self,):
        pass
    
    resultIndices = {'stress': lambda nStateVarsUmat,kw, : np.arange(6) + (kw['location'])*nStateVarsUmat ,
                    'strain': lambda nStateVarsUmat,kw, : np.arange(6) + (kw['location'])*nStateVarsUmat +6,
                    'sdv':    lambda nStateVarsUmat,kw, : kw['indices'] +((kw['location'])-1)*nStateVarsUmat }  
                    
    def getResult(self, **kw):    
        stateVarIndices = self.resultIndices[kw['result']](self.nStateVarsUmat, kw)
        return np.asarray(self.stateVars)[stateVarIndices]
    
    def __dealloc__(self):
        del self.backendElement
