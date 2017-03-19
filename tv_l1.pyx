#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:25:25 2017

@author: andyjones
"""

import numpy as np

cdef extern int tautString_TV1(double *y,double lambd,double *x,int n)
cdef extern void prox_dp(double *y, double lam, double *beta, int n)


def solve(double[:, :] arr, double[:] lambdas):
    cdef int i
    cdef int height = arr.shape[0]
    cdef int width = arr.shape[1]
    cdef double[:, :] output = np.zeros((height, width), dtype=np.float64)
    
    if width > 0:     
        for i in range(height):
            prox_dp(&arr[i, 0], lambdas[i], &output[i, 0], width)
    
    return np.asarray(output)


