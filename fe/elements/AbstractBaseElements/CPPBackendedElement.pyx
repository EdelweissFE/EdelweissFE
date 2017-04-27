#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:35:06 2017

@author: matthias
"""

cdef class BackendedElement:
    
    cdef NISTParallelizableBackendElement* getBackendElement(self,):
        return self.backendElement
    pass
#    cdef NISTParallelizableBackendElement* backendElement