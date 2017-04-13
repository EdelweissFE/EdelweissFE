#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  21 20:17:58 2017

@author: matthias
"""
from libcpp.string cimport string


cdef public bint notificationToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return True
    
cdef public bint warningToMSG(const string cppString):
#    print(cppString.decode('UTF-8'))
    return False

cdef extern from "umatModLeon.h":
    void umatModLeon(       double[],double[],double[],double&,double&,
                            double&,double&,double[],double[],double&, 
                            const double*,const double*,const double*,
                            const double&,const double&,const double&,
                            const double*,const double*,
                            const char*,
                            const int&,const int&,const int&,
                            const int&,const double*,const int&,
                            const double*,const double*,double&,
                            const double&,const double*,const double*,
                            const int&,const int&,const int&,const int&,
                            const int*,const int&,const int)
    
cdef extern from "umatModLeonNonLocal.h":
    void umatModLeonNonLocal(double[],double[],double[],double&,double&,
                            double&,double&,double[],double[],double&, 
                            const double*,const double*,const double*,
                            const double&,const double&,const double&,
                            const double*,const double*,
                            const char*,
                            const int&,const int&,const int&,
                            const int&,const double*,const int&,
                            const double*,const double*,double&,
                            const double&,const double*,const double*,
                            const int&,const int&,const int&,const int&,
                            const int*,const int&,const int)
    
cdef extern from "umatLinearElastic.h":
    void umatLinearElastic(double[],double[],double[],double&,double&,
                            double&,double&,double[],double[],double&, 
                            const double*,const double*,const double*,
                            const double&,const double&,const double&,
                            const double*,const double*,
                            const char*,
                            const int&,const int&,const int&,
                            const int&,const double*,const int&,
                            const double*,const double*,double&,
                            const double&,const double*,const double*,
                            const int&,const int&,const int&,const int&,
                            const int*,const int&,const int)
        
#cdef extern from "umatModLeonPS.h":
#     void umatModLeonPS(     double[],double[],double[],double&,double&,
#                             double&,double&,double[],double[],double&, 
#                             const double*,const double*,const double*,
#                             const double&,const double&,const double&,
#                             const double*,const double*,
#                             const char*,
#                             const int&,const int&,const int&,
#                             const int&,const double*,const int&,
#                             const double*,const double*,double&,
#                             const double&,const double*,const double*,
#                             const int&,const int&,const int&,const int&,
#                             const int*,const int&,const int)    

cdef umatType getUmat(str name):
    
    if name == "modleon":
        return umatModLeon
    if name == "linearElastic":
        return umatLinearElastic
    # elif name == "modleonplanestress":
        # return umatModLeonPS
    elif name == "modleonnonlocal":
        return umatModLeonNonLocal
    
