#/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linspace
from math import pi, sin, factorial
import matplotlib.pyplot as plt
from functools import reduce
import operator

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

def gen_sin_newton(n):
    """Return a function which interpolates sinus"""

    xs = [2*pi*x for x in linspace(0, 1, n)]
    ys = [sin(x) for x in xs]

    (diffs, _) = div_diff(xs, ys)

    return lambda x: horner(float(x), xs, diffs)

def f_derived(n):
    """Return the n-th derivation of sinus.

    Only works for the sinus we do have here, really."""
    def base(n, x):
        return 2**n * pi**n * sin(2*pi*x)

    if (n / 2.0) % 2 == 0:
        return lambda x: -base(n, x)
    else:
        return lambda x: base(n, x)

def cubic_splines(xs, ys):
    if len(xs) != len(ys):
        raise ValueError, "len(xs) != len(ys)"

    n = len(xs) - 1 # degrees of freedom
    _, diffs = div_diff(xs, ys)
    print diffs

    def mu(i):
        return (xs[i] - xs[i-1])/(xs[i+1] - x[i-1])


# Create out own sinus function
functions = [(n, gen_sin_newton(n)) for n in range(5, 20, 2)]
points = linspace(0, 2*pi, 500)

## Calculate the maximum error for each of the interpolations
Es = []
for (n, my_sin) in functions:
    Es.append((n, max([abs(sin(x) - my_sin(x)) for x in points])))

## Plot the results
plt.semilogy(*zip(*Es))

Es = []
## Let's calculate the guesses
for (n, _) in functions:
    xs = linspace(0, 1, n)
    fd = f_derived(n+1)
    errors = []
    for p in points:
        fres = max([fd(e) for e in linspace(0, 1, 500)])
        prod = reduce(operator.mul, [abs(p - xs[i]) for i in range(n)], 1)
        errors.append((fres/factorial(n+1)) * prod)

    Es.append((n, max(errors)))

plt.semilogy(*zip(*Es))

plt.xlabel(u'Anzahl der Interpolationspunkte')
plt.ylabel(u'Maximaler Feheler')
plt.legend((u'Newton', u'Abschätzung'))
plt.title(ur'Interpolationsfehler von $\sin(2\pi x)$')
plt.show()
