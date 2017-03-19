#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:25:25 2017

@author: andyjones
"""

import numpy as np

cdef extern void johnson(double *y, double *beta, int n, double lam)

def solve(double[:, :] arr, double[:] lambdas):
    cdef int i
    cdef int height = arr.shape[0]
    cdef int width = arr.shape[1]
    cdef double[:, :] y = np.zeros((height, width+1), dtype=np.float64)
    cdef double[:, :] beta = np.zeros((height, width+1), dtype=np.float64)
    
    y[:, :-1] = arr
    
    if width > 0:     
        for i in range(height):
            johnson(&y[i, 0], &beta[i, 0], width+1, lambdas[i])
    
    return np.asarray(beta[:, :-1])


