#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 09:16:51 2017

@author: andyjones
"""
import seaborn as sns
import matplotlib.pyplot as plt
import tv_l1
import scipy as sp
import cvxpy as cp

def run(lambd=1):
    signal = sp.zeros((50,))
    signal[25:] = 2
    sp.random.seed(0)
    noise = sp.random.normal(size=50)
    ys = signal + noise
    xs = tv_l1.solve(ys[None, :], sp.array([lambd], dtype=float))[0]
    
    plt.plot(ys, linestyle='none', marker='o')
    plt.plot(signal)
    plt.plot(xs)
    
#run(10)

def quad_solve(ys, lambd):
    xs = cp.Variable(len(ys))
    
    loss = 0.5*cp.sum_entries(cp.square(ys - xs)) + lambd*cp.norm1(cp.diff(xs)) + lambd*cp.norm1(xs[-1])
    cp.Problem(cp.Minimize(loss)).solve()
    
    return xs.value.A[:, 0]

def example_tube(lambd=1):
    sp.random.seed(0)
    ys = sp.array([-1, -3, -2], dtype=float)
    n = len(ys)
    rs = sp.cumsum(ys)
    
    tube_mid = sp.hstack([[0], rs])
    tube_lower = tube_mid.copy()
    tube_lower[1:-1] -= lambd
    tube_upper = tube_mid.copy()
    tube_upper[1:-1] += lambd
    
    fig, (top, bot) = plt.subplots(2, sharex=True)
    top.plot(tube_mid)
    top.fill_between(sp.arange(n+1), tube_lower, tube_upper, alpha=0.3)
    top.set_xlim(0, n)
    
    xs = tv_l1.solve(ys[None, :].copy(), sp.array([lambd], dtype=float))[0]
    top.plot(sp.hstack([[0], sp.cumsum(xs)]))
    
    bot.plot(sp.hstack([[sp.nan], ys]))
    bot.plot(sp.hstack([[sp.nan], xs]))
    
    exact = quad_solve(ys, lambd)
    top.plot(sp.hstack([[0], sp.cumsum(exact)]))
    
    bot.plot(sp.hstack([[sp.nan], exact]))
    
    top.legend(['input', 'fast', 'quad'], loc='top-left')
    
    
example_tube(5)
    