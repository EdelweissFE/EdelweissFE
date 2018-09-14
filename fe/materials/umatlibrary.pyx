#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  21 20:17:58 2017

@author: matthias
"""
from libcpp.string cimport string

cdef public bint notificationToMSG(const string* cppString):
    # printing is not possible, as the GIL is usually released for parallel computing
    return True
    
cdef public bint warningToMSG(const string cppString):
    # printing is not possible, as the GIL is usually released for parallel computing
    return False

