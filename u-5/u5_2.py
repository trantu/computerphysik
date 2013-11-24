#/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linspace
from math import pi, sin, cos, factorial
import matplotlib.pyplot as plt
from functools import reduce
import operator
from copy import deepcopy

def div_diff(xs, starting_values):
    """Calculate the Newton dividing differences.

    Return the "top layer" which are the f(x0,...,xi) values
    """

    n = len(xs)
    fprev = starting_values
    fhist = []
    ret = []

    for k in range(1, n+1):
        fhist.append(fprev)
        ret.append(fprev[0])
        fprev = [(fprev[i+1] - fprev[i]) / (xs[k+i] - xs[i]) for i in range(n-k)]

    return (ret, fhist)

def horner(x, xs, diffs):
    n = len(xs)
    s = diffs[n-1]
    for k in range(n-2, -1, -1):
        s = s * (x - xs[k]) + diffs[k]

    return s

def newton_inter(xs, ys, x):
    """Interpolate f such that ys = f(xs) at the point x. """

    if len(xs) != len(ys):
        raise ValueError, "len(xs) != len(ys)"

    # Let's start with the dividing differences
    (diffs, _) = div_diff(xs, ys)
    
    # and use the Horner schema to interpolate
    return horner(x, xs, diffs)

class NewtonInterpolation:
    """A function which interpolates sinus with 'n' interpolation points"""
    def __init__(self, xs, ys):
        if len(xs) != len(ys):
            raise ValueError, "len(xs) != len(ys)"

        self.xs = xs
        (self.diffs, _) = div_diff(xs, ys)

    def __call__(self, x):
        return horner(float(x), self.xs, self.diffs)

def f_derived(n):
    """Return the n-th derivation of sinus.

    Only works for the sinus we do have here, really."""
    def base(n, x):
        return 2**n * pi**n * cos(2*pi*x)

    if (n / 2) % 2 == 0:
        return lambda x: -base(n, x)
    else:
        return lambda x: base(n, x)


class CubicSpline:
    def __init__(self, xs, ys):
        if len(xs) != len(ys):
            raise ValueError, "len(xs) != len(ys)"

        self.xs = xs
        self.ys = ys
        self.moments()

    def moments(self):
        xs = self.xs
        ys = self.ys

        n = len(xs) - 1 # degrees of freedom

        mu = [0] * (n+1)
        for i in range(1, n):
            mu[i] = (xs[i] - xs[i-1])/(xs[i+1] - xs[i-1])
            mu[n] = (xs[n] - xs[n-1])/(xs[n] - xs[1] - xs[0] - xs[n-1])

        A = np.matrix(np.zeros((n-1, n-1)))
        #print "n", n, "A", A
        for i in range(n-1):
            if i != 0:
                A[i,i-1] = mu[i+1]
            A[i,i] = 2
            if i != n-2:
                A[i,i+1] = 1.0 - mu[i+1]

        b = np.matrix(np.zeros((n-1, 1)))
        (d, dd) = div_diff(xs, ys)
        #print "d", d, "dd", dd
        for i in range(0, n-1):
            b[i,0] = 6.0*dd[2][i]

        [self.M] = np.linalg.solve(A, b).flatten().tolist()

    def __call__(self, x):
        xs = self.xs
        ys = self.ys
        n = len(xs)
        # Add x_{n+1} and y)_{n+1} to get the edge cases
        #xs2 = xs + [xs[-1] + xs[1] + xs[0]]
        #ys2 = ys + [ys[-1]]

        M = [0.0] + self.M + [0.0]
        #M.insert(0, 0.0)
        #M.append(0.0)

        h = [xs[i] - xs[i-1] for i in range(n)]

        C = [0] * (n)
        for i in range(n):
            first = (ys[i] - ys[i-1]) / h[i]
            second = (h[i]/6.0) * (M[i] - M[i-1])
            C[i] = first - second

        # def C(i):
            # first = (ys[i] - ys[i-1]) / h[i]
            # second = (h[i]/6.0) * (M[i] - M[i-1])
            # return first - second

        D = [0] * n
        for i in range(n):
            first = (ys[i] + ys[i-1]) / 2.0
            second = ((h[i]**2)/12.0) * (M[i] + M[i-1])
            D[i] = first - second

        # def D(i):
        #     first = (ys[i] + ys[i-1]) / 2.0
        #     second = ((h[i]**2)/12.0) * (M[i] + M[i-1])
        #     return first - second

        def s(x, i):
            first = M[i-1] * (((xs[i] - x)**3) / (6*h[i]))
            second = M[i] * (((x - xs[i-1])**3) / (6*h[i]))
            third = C[i] * (x - ((xs[i-1]+xs[i])/2.0))
            fourth = D[i]
        #print 'first', first, 'second', second, 'third', third, 'fourth', fourth
            return first + second + third + fourth

        # find what interval x is in
        i = -1
        for ii in range(1, n):
            if x >= xs[ii-1] and x <= xs[ii]:
                i = ii
                break

        if i == -1:
            raise ValueError, ("No such interval for %f" % x)

        return s(x, i)
        

ns = range(5, 20, 2)
## Calculate the maximum error for each of the interpolations
Es = []
for n in ns:
    xs = linspace(0, 1, n)
    ys = [sin(2*pi*x) for x in xs]
    my_sin = NewtonInterpolation(xs, ys) # our sin() for this round
    error = max([abs(sin(2*pi*x) - my_sin(x)) for x in linspace(0, 1, 500)])
    Es.append(error)

## Plot the results
plt.semilogy(ns, Es)

# Let's do cublic splines
Es = []
for n in ns:
    xs = linspace(0, 1, n)
    print "xs", xs
    ys = [sin(2*pi*x) for x in xs]
    my_sin = CubicSpline(xs, ys)
    error = max([abs(my_sin(p) - sin(2*pi*p)) for p in linspace(0, 1, 500)])
    Es.append(error)

plt.semilogy(ns, Es)

# 5.2.3 Bestimmung des Maximalfehlers durch $\frac{f^{n+1}(\xi)}{(n+1)!} \prod^n_{i=0}(x-x_i)$.
#
# Wir befinden uns im Intervall [0,1], also sind x wie x_i immer <= 1,
# womit (x - x_i) auch immer <= 1 ist und genau so ist der Produkt <=
# 1. Das können wir also aus der Formel weglassen. Die (n+1)-te
# Ableitung von sin() (mit (n+1) eine ungerade Zahl, wie es in diesem
# Fall immer ist) ist abs(2**(n+1) * pi**(n+1) * cos(2*pi*x)). Der
# Cosinus ist wiederum maximal eins und wir können ihn ignorieren.
#
# Es bleibt also 2**(n+1) * pi**(n+1) / (n + 1)! als maximaler Fehler
# für ein bestimmtes n. Hier haben wir n als Grad des Polynoms
# benutzt, im Programm ist aber n die Anzahl der Stützpunkte. Also
# wird im Code n geschrieben, statt (n+1) wie auf der mathematische
# Formel.

Es = [(2**n * pi**n) / factorial(n) for n in ns]

plt.semilogy(ns, Es)

plt.xlabel(u'Anzahl der Interpolationspunkte')
plt.ylabel(u'Maximaler Feheler')
plt.legend((u'Newton', u'splines', u'Abschätzung'))
plt.title(ur'Interpolationsfehler von $\sin(2\pi x)$')
plt.show()
