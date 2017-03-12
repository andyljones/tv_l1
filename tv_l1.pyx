#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:25:25 2017

@author: andyjones
"""

import numpy as np

cdef extern void TV1D_denoise(double* arr, double* output, const int width, const double lambd)

def solve(double[:, :] arr, double[:] lambdas):
    cdef int i
    cdef int height = arr.shape[0]
    cdef int width = arr.shape[1]
    cdef double[:, :] output = np.zeros((height, width), dtype=np.float64)
    
    if width > 0:     
        for i in range(height):
            TV1D_denoise(&arr[i, 0], &output[i, 0], width, lambdas[i])
    
    return np.asarray(output)


