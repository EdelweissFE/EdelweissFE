#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 08:50:46 2017

@author: matthias
"""

import numpy as np

def sigPrinc(x):
    if x.ndim==2:
        x_ = x
        eig = np.asarray([ np.linalg.eig( np.array([[x[0], x[3], x[4]],
                                                  [x[3], x[1], x[5]],
                                                  [x[4], x[5], x[2]]])  )[0] for x in x_])
        return eig

sympyMathModules = {'sigPrinc':sigPrinc}