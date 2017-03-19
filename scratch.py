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
    
    loss = 0.5*cp.sum_entries(cp.square(y - x)) + lambd*cp.norm1(cp.diff(x))
    cp.Problem(cp.Minimize(loss), [x[0] == y[0]]).solve()
    
    return x.value.A[:, 0]
    
def generate_example():
    n = sp.random.randint(3, 20)
    return {'y': sp.random.normal(scale=1, size=(n,)),
            'lambd': sp.random.uniform(.1, 2)}

def plot_example(y, lambd):
    quad = quad_solve(y, lambd)
    taut = tv_l1.solve(y[None, :].astype(float), sp.array([lambd], dtype=float))[0]
    
    n = len(y)
    fig, bot = plt.subplots(1, sharex=True)    
    bot.set_xlim(0, n)
    bot.set_ylim(y.min()-1, y.max()+1)
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], y]), label='input')
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], quad]), label='quad opt soln', marker='o')
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], taut]), label='johnson soln')
    bot.legend()
    
#example = generate_example()
#plot_example(**example)