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
    
    loss = 0.5*cp.sum_entries(cp.square(y - x)) + lambd*cp.norm1(cp.diff(x))#ß + lambd*cp.norm1(x[-1])
    cp.Problem(cp.Minimize(loss)).solve()
    
    return x.value.A[:, 0]

def tautstring(y, lambd):
    y = y.copy()
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
        lambd_hat = lambd if main else 0

        debug.append({k: v for k, v in locals().items() if sp.isscalar(v)})
        
        mn_height += mn - y[i]
        mx_height += mx - y[i]

        if mn_height > +lambd_hat:
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
            
def dyn_prog(y, lambd):
    N = len(y)
    if N == 0:
        return sp.array([])
    if (N == 1) or (lambd == 0):
        return y.copy()
    
    beta = sp.nan*sp.zeros(N)
        
    x = sp.nan*sp.zeros(2*N)
    a = sp.nan*sp.zeros(2*N)
    b = sp.nan*sp.zeros(2*N)
    
    t_m, t_p = sp.nan*sp.zeros(N-1), sp.nan*sp.zeros(N-1)
    
    # First iteration
    a_lo, a_hi = +1, -1
    b_lo, b_hi = -y[0], +y[0] 
    t_m[0], t_p[0] = y[0]-lambd, y[0]+lambd
    l, r = N-1, N
    x[l], x[r] = t_m[0], t_p[0]

    a[l], a[r] = a_lo, a_hi
    b[l], b[r] = b_lo+lambd, b_hi+lambd

    a_first, a_last = +1, -1
    b_first, b_last = -y[1]-lambd, +y[1]-lambd

    debug = []

    # Main loop
    for k in range(1, N-1):
        debug.append({k: v for k, v in locals().items() if sp.isscalar(v)})
        
        # Compute negative knot
        a_lo, b_lo = a_first, b_first
        for lo in range(l, r+1):
            print('lo', lo)
            if (+a_lo*x[lo]+b_lo > -lambd):
                l = lo-1
                break
            a_lo += a[lo]
            b_lo += b[lo]
        else:
            l = r
        t_m[k] = (-b_lo-lambd)/(+a_lo)
        print('l', l)
        x[l] = t_m[k]

        # Compute positive knot
        a_hi, b_hi = a_last, b_last
        for hi in range(r, l-1, -1):
            print('hi', hi)
            if (-a_hi*x[hi]-b_hi < +lambd):
                r = hi+1
                break
            a_hi += a[hi]
            b_hi += b[hi]
        else:
            r = hi
        t_p[k] = (+b_hi+lambd)/(-a_hi)
        print('r', r)
        x[r] = t_p[k]

        a[l], a[r] = a_lo, a_hi
        b[l], b[r] = b_lo+lambd, b_hi+lambd
        
        a_first, a_last = +1, -1
        b_first, b_last = -y[k+1]-lambd, +y[k+1]-lambd

    debug.append({k: v for k, v in locals().items() if sp.isscalar(v)})
        

    # Last coefficient
    a_lo, b_lo = a_first, b_first
    for lo in range(l, r+1):
        if (+a_lo*x[lo]+b_lo > 0):
            break
        a_lo += a[lo]
        b_lo += b[lo]
    beta[N-1] = -b_lo/(+a_lo)

    # Rest of the coeffs
    for k in range(N-2, -1, -1):
        if beta[k+1] > t_p[k]:
            beta[k] = t_p[k]
        elif beta[k+1] < t_m[k]:
            beta[k] = t_m[k]
        else:
            beta[k] = beta[k+1]

    return beta, pd.DataFrame(debug)

        
    
def generate_example():
    n = sp.random.randint(3, 20)
    return {'y': sp.random.normal(scale=1, size=(n,)),
            'lambd': sp.random.uniform(.1, 5)}
    
def test_example(y, lambd):
    quad = quad_solve(y, lambd)
    taut, debug = dyn_prog(y, lambd)
    
    return sp.allclose(quad, taut), quad, taut

def plot_example(y, lambd):
    quad = quad_solve(y, lambd)
    taut, debug = dyn_prog(y, lambd)
    print(debug)
#    taut = tv_l1.solve(y[None, :].astype(float), sp.array([lambd], dtype=float))[0]
    
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

plot_example(**generate_example())
#y = sp.array([-2, 2, 2])
#lambd = 1
#plot_example(y, lambd)