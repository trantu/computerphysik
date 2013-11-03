#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

import math
from scipy.misc import derivative
import matplotlib.pyplot as plot
from operator import sub

# return 1 if v is positive, 0 otherwise
def sig(v):
    return v >= 0

def next_interval(f, (a,b)):
    """
    Return a more accurate interval in which the root can be found
    """

    fa = f(a)
    fb = f(b)

    if fa * fb >= 0:
        raise ValueError("there is no root in this interval")

    mid = (a+b) / 2
    fmid = f(mid)

    if sig(fa) != sig(fmid):
        return (a, mid)
    else:
        return (mid, b)

def bisection(f, a, b, delta):
    """
    Find a root for `f` in the interval (a, b)
    """

    # Start with the user's interval
    interval = (a, b)

    # Keep getting more accurate invervals until the differences are
    # smaller than delta
    n = 0
    while not abs(sub(*interval)) < delta:
        interval = next_interval(f, interval)
        n += 1

    return (interval, n)

def fixed_point_iteration(f, x0, delta):
    """
    Return the zero with `delta` precision
    """

    last_point = x0
    point = f(x0)

    n = 0
    while not abs(point - last_point) < delta:
        last_point = point
        point = f(last_point)
        n += 1

    return (point, n)

def next_point(f, x):
    """ Return the next point for Newton """
    deriv = derivative(f, x)
    return x - (f(x)/deriv)

def newton(f, x0, delta):
    """
    Return the zero with `delta` precision
    """

    last_point = x0
    point = next_point(f, x0)

    n = 0
    while not abs(last_point - point) < delta:
        last_point = point
        point = next_point(f, last_point)
        n += 1

    return (point, n)

## Sanity checking
f = lambda x: math.cos(x)**2 - x + 0.2
fi = lambda x: math.cos(x)**2 + 0.2
# print("For f(x):")
# ((a,b), n) = bisection(f, -1.0, 1.0, 0.00001)
# print("\tbisection (%g, %g), %g" % (a,b,n))
# print("\tfixed point iter %g, %d" % fixed_point_iteration(fi, 0.0, 0.01))
# print("\tNewton iter %g, %d" % newton(f, 0.0, 0.01))

g = lambda x: x**2 - x - 0.2
gi = lambda x: x**2 - 0.2
# print("For g(x):")
# #print("\tbisection (%g, %g)" % bisection(g, -1.0, 1.0, 0.01))
# print("\tfixed point iter %g, %d" % fixed_point_iteration(gi, 0.0, 0.0000000001))
# print("\tNewton iter %g, %d" % newton(g, 0.0, 0.000000000000001))

deltas = [2**-x for x in range(1, 21)]

for (fn, fni, t) in [(f, fi, u'$f(x)$'), (g, gi, u'$f(x)$')]:
    steps_bisection = [n for (_, n) in map(lambda d: bisection(fn, -1.0, 1.0, d), deltas)]
    steps_fixed_point = [n for (_, n) in
                         map(lambda d: fixed_point_iteration(fni, 0.0, d), deltas)]
    steps_newton = [n for (_, n) in map(lambda d: newton(fn, 0.0, d), deltas)]

    plot.title(u'n gegen $\Delta$ bei %s' % t)
    plot.loglog(deltas, steps_bisection)
    plot.loglog(deltas, steps_fixed_point)
    plot.loglog(deltas, steps_newton)
   
    plot.ylabel(u'n')
    plot.xlabel(u'$\Delta$')
    plot.legend((u'Bisektion', u'Fixpunkt', u'Newton'))

    plot.show()
