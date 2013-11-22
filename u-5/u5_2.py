#/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linspace
from math import pi, sin
import matplotlib.pyplot as plt
from functools import reduce

def div_diff(xs, starting_values):
    """Calculate the Newton dividing differences.

    Return the "top layer" which are the f(x0,...,xi) values
    """

    n = len(xs)
    fprev = starting_values
    ret = []

    for k in range(1, n+1):
        ret.append(fprev[0])
        fs = []
        fprev = [(fprev[i+1] - fprev[i]) / (xs[k+i] - xs[i]) for i in range(n-k)]

    return ret

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
    diffs = div_diff(xs, ys)
    
    # and use the Horner schema to interpolate
    return horner(x, xs, diffs)

def gen_sin_newton(n):
    """Return a function which interpolates sinus"""

    xs = [2*pi*x for x in linspace(0, 1, n)]
    ys = [sin(x) for x in xs]

    diffs = div_diff(xs, ys)

    return lambda x: horner(float(x), xs, diffs)


# Create out own sinus function
functions = [(n, gen_sin_newton(n)) for n in range(5, 20, 2)]
points = linspace(0, 4*pi, 500)

## Calculate the maximum error for each of the interpolations
Es = []
for (n, my_sin) in functions:
    Es.append((n, max([abs(sin(x) - my_sin(x)) for x in points])))

## Plot the results
plt.semilogy(*zip(*Es))
plt.xlabel(u'Interpolationspunkte')
plt.ylabel(u'Maximaler Feheler')
plt.legend((u'Newton', ))
plt.title(ur'Interpolation von $\sin(2\pi x)$')
plt.show()
