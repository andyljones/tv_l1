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

def quad_solve(y, lambd):
    x = cp.Variable(len(y))
    
    loss = 0.5*cp.sum_entries(cp.square(y - x)) + lambd*cp.norm1(cp.diff(x)) + lambd*cp.norm1(x[-1])
    cp.Problem(cp.Minimize(loss)).solve()
    
    return x.value.A[:, 0]

def tautstring(y, lambd):
    N = len(y)
    
    mn_height, mx_height = 0, 0
    mn, mx = y[0] - lambd, y[0] + lambd
    
    i = 0
    last_break = -1
    mn_break, mx_break = 0, 0
    
    x = sp.nan*sp.zeros_like(y)
    
    debug = []
#    import pdb; pdb.set_trace()
    while i < N:
        main = (i < N-1)
        lambd_hat = lambd if main else lambd

        debug.append({k: v for k, v in locals().items() if sp.isscalar(v)})
        
        mn_height += mn - y[i]
        mx_height += mx - y[i]

        if  mn_height > +lambd_hat:
            i = mn_break + 1
            x[last_break+1:mn_break+1] = mn
            last_break = mn_break
            mn, mx = y[i], y[i] + 2*lambd
            mn_height, mx_height = (-lambd, +lambd) if main else (-lambd, -lambd)
            mn_break, mx_break = i, i
            if main:
                i += 1    
        elif mx_height < -lambd_hat:
            i = mx_break + 1
            x[last_break+1:mx_break+1] = mx
            last_break = mx_break
            mn, mx = y[i] - 2*lambd, y[i]
            mn_height, mx_height = (-lambd, +lambd) if main else (+lambd, +lambd)
            mn_break, mx_break = i, i
            if main:
                i += 1
        else:
            if mx_height > +lambd_hat:
                mx += (+lambd_hat - mx_height)/(i - last_break)
                mx_height = +lambd_hat
                mx_break = i
                
            if mn_height < -lambd_hat:
                mn += (-lambd_hat - mn_height)/(i - last_break)
                mn_height = -lambd_hat
                mn_break = i
            
            i += 1
    
    
    debug.append({k: v for k, v in locals().items() if sp.isscalar(v)})
    
    
    x[last_break+1:] = mn if abs(mn) < abs(mx) else mx

    return x, pd.DataFrame(debug)
            
def generate_example():
    n = sp.random.randint(2, 50)
    return {'y': sp.random.normal(size=(n,)),
            'lambd': sp.random.uniform(.1, 1)}
    
def test_example(y, lambd):
    quad = quad_solve(y, lambd)
    taut, debug = tautstring(y, lambd)
    
    return sp.allclose(quad, taut), quad, taut

def plot_example(y, lambd):
    quad = quad_solve(y, lambd)
    taut, debug = tautstring(y, lambd)
    
    n = len(y)
    r = sp.cumsum(y)
    
    tube_mid = sp.hstack([[0], r])
    tube_lower = tube_mid.copy()
    tube_lower[1:-1] -= lambd
    tube_upper = tube_mid.copy()
    tube_upper[1:-1] += lambd

    fig, (top, bot) = plt.subplots(2, sharex=True)

    top.fill_between(sp.arange(n+1), tube_lower, tube_upper, alpha=0.3)    
    top.plot(tube_mid, label='input')
    top.plot(sp.hstack([[0], sp.cumsum(quad)]), label='quad', marker='o')
    top.plot(sp.hstack([[0], sp.cumsum(taut)]), label='taut')
    top.set_xlim(0, n)
    top.legend(loc='top left')
    
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], y]), label='input')
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], quad]), label='quad', marker='o')
    bot.step(sp.arange(n+1), sp.hstack([[sp.nan], taut]), label='taut')

plot_example(**generate_example())