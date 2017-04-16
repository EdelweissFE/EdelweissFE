#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 13:54:09 2017

@author: matthias
"""

from fe.materials.umatlibrary cimport pUmatType 

ctypedef void (*pSimpleUelWithUmatType) ( double[], double[], double[], const int &, const double*, 
              const int &, const double*, const double*, const double*, 
              const double*, const double &, const int& , double &, 
              const int*, const int&, pUmatType, int ) nogil