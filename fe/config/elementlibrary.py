#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  17 19:10:42 2017

@author: matthias
"""
from  fe.elements.uelcps4.element import Element as uelCPS4
from  fe.elements.uelcps4nonlocal.element import Element as uelCPS4NonLocal

elementlibrary = {"uelCPS4" :         uelCPS4,
                  "uelCPS4NonLocal" : uelCPS4NonLocal,}