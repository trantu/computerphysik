#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import exp
import numpy as np
import scipy.integrate

from itertools import islice

def take(seq, n):
    return list(islice(seq, 0, n))

def V(x):
    return x**8 * 2.5e-4

def Valt(x):
    return 0.5 * x**2

def schroed_system(E, V):
    def func(y, x):
        a, b = y
        return [b, -2.0 * a * (E - V(x))]

    return func

def gradient(E, V):
    """Ignore this"""
    def func(y, t):
        #return [[-2*(E - V(x)), 0], [0, 1]]
        return [[0, 1], [-2*(E - V(x)), 0]]

def one_step(func, y0, rng):
    """Return the next value we're interested in for the Eigenvalues"""
    return scipy.integrate.odeint(func, y0, rng)[-1]

def sign_change(a, b):
    return True if a * b < 0 else False

def find_eigenvalues(E0, V, y0, rng, step=0.01, delta=0.001):
    # Set up the initial values
    func = schroed_system(E0, V)
    [E, _] = one_step(func, y0, rng)
    if E < delta:
        yield E

    # We're looking for places where \Phi(x)=0 so we keep searching
    # for a change in the sign
    Eprev = E
    while True:
        E = one_step(func, y0, rng)
        if sign_change(Eprev, E):
            # bisect a more exact number
            # yield that number

rng = np.arange(-6, 6, 0.01)
#y0 = [10e-10, 10e-10]
y0 = [10e-7, 10e-7]
EV = (0.5001, Valt)

find_eigenvalues(schroed_system(*EV), y0, rng)
