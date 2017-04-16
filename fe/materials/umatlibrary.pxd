#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  21 20:28:33 2017

@author: matthias
"""

ctypedef void (*pUmatType)( double[],double[],double[],double&,double&,
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

cdef pUmatType getUmat(str name)
