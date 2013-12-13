#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import exp
import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt

def euler(n, h, y0, F, t0=0):
    y, tn = float(y0), float(t0)
    ys = []
    print 'h', h
    for i in range(1, n+1):
        y = y + h*F(tn, y)
        #print 'euler y', y, 'tn', tn
        ys.append(y)
        tn = tn + h

    return ys

def mittelpunktsregel(n, h, y0, F, t0=0):
    y, tn = float(y0), float(t0)
    hh = h/2.0
    ys = []
    for i in range(1, n+1):
        y = y + h*F(tn + hh, y + hh*f(tn, y))
        ys.append(y)
        tn = tn + h

    return ys

def runge_kutta(n, h, y0, F, t0=0):
    y, tn = float(y0), float(t0)
    hh = h/2.0
    ys = []
    for i in range(1, n+1):
        K1 = F(tn, y)
        K2 = F(tn + hh, y + hh*K1)
        K3 = F(tn + hh, y + hh*K2)
        K4 = F(tn + h,  y + h*K3)
        y = y + (h/6.0)*(K1 + 2*K2 + 2*K3 + K4)
        ys.append(y)
        print 'tn', tn, 'y', ys[-1]
        tn = tn + h

    return ys

def exact(n, h, y0, f, t0=0):
    y, tn = float(y0), float(t0)
    ys = []
    for i in range(1, n+1):
        ys.append(np.exp(-tn))
        tn = tn + h
    return ys

## Trivial test
#assert(16 == euler(4, 1, 1, lambda t, y: y)[-1])

trange = (0, 2)
hs = [0.1, 0.01, 0.001]
funs = [
    ('Euler', euler, '-'),
    ('Mittelpunkt', mittelpunktsregel, '-'),
    ('Runge-Kutta', runge_kutta, '-'),
    #('Exakt', exact, '--')
]

def f(t, y):
    return 3.0*y - 4.0*exp(-t)

for h in hs:
    for (lbl, fn, opt) in funs:
        # calcultate n for the range
        n = int((trange[1] - trange[0])/h) + 1
        ts = linspace(0, 2, n)
        ys = fn(n, h, 1, f)
        plt.plot(ts, ys, opt, label=lbl + ' h = %f' % h)
    plt.legend()

# calcultate n for the range
h = 0.1
n = int((trange[1] - trange[0])/h) + 1
ts = linspace(0, 2, n)
ys = euler(n, h, 1, f)
plt.plot(ts, ys, '--', label='Exakt')

plt.show()
        
    
