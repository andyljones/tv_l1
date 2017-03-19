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
import pandas as pd

def quad_solve(y, lambd):
    x = cp.Variable(len(y))
    
    loss = 0.5*cp.sum_entries(cp.square(y - x)) + lambd*cp.norm1(cp.diff(x)) + lambd*cp.norm1(x[-1])
    cp.Problem(cp.Minimize(loss)).solve()
    
    return x.value.A[:, 0]
    
def generate_example():
    n = sp.random.randint(3, 20)
    return {'y': sp.random.normal(scale=1, size=(n,)),
            'lambd': sp.random.uniform(.1, 2)}

def plot_example(y, lambd):
    quad = quad_solve(y, lambd)
    taut = tv_l1.solve(y[None, :].astype(float), sp.array([lambd], dtype=float))[0]
    
    n = len(y)
    r = sp.cumsum(y)
    
    tube_mid = sp.hstack([[0], r])
    tube_lower = tube_mid.copy()
    tube_lower[1:] -= lambd
    tube_upper = tube_mid.copy()
    tube_upper[1:] += lambd

    fig, (top, bot) = plt.subplots(2, sharex=True)

    top.fill_between(sp.arange(n+1), tube_lower, tube_upper, alpha=0.3)    
    top.plot(tube_mid, label='input')
    top.plot(sp.hstack([[0], sp.cumsum(quad)]), label='quad', marker='o')
    top.plot(sp.hstack([[0], sp.cumsum(taut)]), label='taut')
    top.set_xlim(0, n)
    top.legend(loc='upper left')
    
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], y]), label='input')
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], quad]), label='quad', marker='o')
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], taut]), label='taut')

#plot_example(**generate_example())
y = sp.array([-2, 2, 2])
lambd = 1
plot_example(y, lambd)