#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:25:25 2017

@author: andyjones
"""

import numpy as np
cimport cython

cdef extern void johnson(double *y, double *beta, int n, double lam)

@cython.boundscheck(False)
def solve(double[:, :] y, double[:] lambdas):
    """For each row $y_i$ and $\lambd_i$, finds the $\beta_i$ which minimizes
    
    $$\min_\beta \frac{1}{2}\|y_i - \beta_i\|^2_2 + \lambda_i | D\beta_i |_1 \text{ s.t } \beta_{i0} = y_{i0}$$
    
    where $D$ is the difference operator. 
    
    Internally this uses a lightly-modified version of Johnson's 2012 dynamic programming algorithm. Johnson's 
    algorithm was more easily modified to solve the above problem when the constraint is that $\beta_{in} = y_{in}$,
    so internally this function reverses each row of $y$, passes it to the algorithm, then reverses the result.
    """
    
    cdef int i
    cdef int n = y.shape[1]
    cdef double[:, :] y_rev = y[:, ::-1].copy()
    cdef double[:, :] beta = y.copy()
    
    for i in range(y.shape[0]):
        johnson(&y_rev[i, 0], &beta[i, 0], n, lambdas[i])
    
    return np.asarray(beta)[:, ::-1]


