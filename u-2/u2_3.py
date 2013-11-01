#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

import math

# return 1 if v is positive, 0 otherwise
def sig(v):
    return v >= 0

def next_interval(f, (a,b)):
    """

    Return a more exact interval in which the zero lays, or the input
    if they're less than `delta` apart.

    """

    fa = f(a)
    fb = f(b)

    if fa * fb >= 0:
        raise ValueError("there is no zero in this interval")

    mid = (a+b) / 2
    fmid = f(mid)

    if sig(fa) != sig(fmid):
        return (a, mid)
    else:
        return (mid, b)

def tuple_diff((a,b)):
    return abs(a - b)

def bisection(f, a, b, delta):
    """
    Find a zero for `f` in the interval (a, b)
    """

    # Start with the user's interval
    interval = (a, b)

    # Keep getting more accurate invervals until the differences are
    # smaller than delta
    while not tuple_diff(interval) < delta:
        interval = next_interval(f, interval)

    return interval

def fixed_point_iteration(f, x0, delta):
    """
    Return the zero with `delta` precision
    """

    last_point = x0
    point = f(x0)

    while not abs(point - last_point) < delta:
        last_point = point
        point = f(last_point)

    return point

f = lambda x: math.cos(x)**2 - x + 0.2
fi = lambda x: math.cos(x)**2 + 0.2
print("bisection (%g, %g)" % bisection(f, -1.0, 1.0, 0.01))
print("fixed point iter %g" % fixed_point_iteration(fi, 0.0, 0.01))
