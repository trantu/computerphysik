#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import exp, sin, pi
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

from itertools import islice, imap

def nth(seq, n):
    return seq[n]

def first(seq):
    return nth(seq, 0)
def second(seq):
    return nth(seq, 1)

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

def schroedinger(E, V, y0, rng):
    func = schroed_system(E, V)
    return scipy.integrate.odeint(func, y0, rng)

def sign_change(a, b):
    return True if a * b < 0 else False

def bisect_root(V, y0, rng, Eleft, E0left, Eright, E0right, delta):
    """Bisect until we find a value smaller than delta"""

    E0mid = (E0left + E0right) / 2.0
    func = schroed_system(E0mid, V)
    [Emid, _] = one_step(func, y0, rng)

    if abs(Emid) < delta:
        return (E0mid, Emid)

    if sign_change(Eleft, Emid):
        return bisect_root(V, y0, rng, Eleft, E0left, Emid, E0mid, delta)
    elif sign_change(Emid, Eright):
        return bisect_root(V, y0, rng, Emid, E0mid, Eright, E0right, delta)
    else:
        assert ValueError, "there is no sign change in the values given"

def find_eigenvalues(E0, V, y0, rng, step=0.01, delta=0.001):
    # Set up the initial values
    func = schroed_system(E0, V)
    [E, _] = one_step(func, y0, rng)
    if E < delta:
        yield (E0, E)

    # We're looking for places where \Phi(x)=0 so we keep searching
    # for a change in the sign and then we bisect to try to find a more
    # exact value
    while True:
        Eprev = E
        E0prev = E0
        E0 = E0 + step
        func = schroed_system(E0, V)
        [E, _] = one_step(func, y0, rng)
        if sign_change(Eprev, E):
            yield bisect_root(V, y0, rng, Eprev, E0prev, E, E0, delta)

rng = np.arange(-6, 6, 0.01)
y0 = [10e-10, 10e-10]
E0 = 0

eigens_gen = find_eigenvalues(E0, Valt, y0, rng)
eigens = take(imap(first, eigens_gen), 5)
print u'Erste fünf Eigenwerte:', eigens

# Now that we have the Eigenvalues, we solve the Scrhödinger equation
# with them.
for E in eigens:
   plt.plot(rng, map(first, schroedinger(E, V, y0, rng)), label='E = %f' % E)

plt.legend()
plt.show()

# 9.2.2

squares = []
for kmax in [10, 40, 100]:
    def Vsquare(x):
        f  = 1.0 / (8.0 / 2.0)
        sums = [sin(2*pi*(2*k - 1)*f*x)/(2*k -1) for k in range(1, kmax+1)]
        return (4.0/pi) * sum(sums)

    plt.figure()
    plt.title('rechteckiger Potenzialtopf $k_{max} = %d$' % kmax)
    vals = [Vsquare(x) for x in rng]
    plt.plot(rng, [Vsquare(x) for x in rng], label='Potenzial')
    squares.append(vals)
    y0 = [1e-10, 1e-10]
    eigens = take(imap(first, find_eigenvalues(E0, Vsquare, y0, rng, 0.01)), 5)
    #print eigens
    for E in eigens:
        plt.plot(rng, map(first, schroedinger(E, Vsquare, y0, rng)), label='E = %f' % E)
    plt.legend()
    
plt.show()
