#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

import math

# return 1 if v is positive, 0 otherwise
def sig(v):
    return v >= 0

def next_interval(f, (a,b), delta):
    """

    Return a more exact interval in which the zero lays, or the input
    if they're less than `delta` apart.

    """

    fa = f(a)
    fb = f(b)

    if abs(fa - fb) < delta:
        return (a, b)

    if fa * fb >= 0:
        raise ValueError("there is no zero in this interval")

    mid = (a+b) / 2
    fmid = f(mid)

    if sig(fa) != sig(fmid):
        return (a, mid)
    else:
        return (mid, b)


def bisection(f, a, b, delta):
    """
    Find a zero for `f` in the interval (a, b)
    """

    # Start with the user's interval
    last_interval = (a, b)
    interval = next_interval(f, last_interval, delta)

    # Keep asking for the next interval as long as we find a
    # change. If they're less than delta apart, next_interval() will
    # return the input
    while interval != last_interval:
        last_interval = interval
        interval = next_interval(f, last_interval, delta)

    return interval

f = lambda x: math.cos(x)**2 - x + 0.2
print("bisection (%g, %g)" % bisection(f, -1.0, 1.0, 0.01))
