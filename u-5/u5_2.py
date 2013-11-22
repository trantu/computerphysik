#/usr/bin/env python
# -*- coding: utf-8 -*-

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

def newton_inter(xs, ys, x):
    """Interpolate f such that ys = f(xs) at the point x. """

    if len(xs) != len(ys):
        raise ValueError, "len(xs) != len(ys)"

    # Let's start with the dividing differences
    n = len(xs)
    diffs = div_diff(xs, ys)
    
    # and use the Horner schema to interpolate
    s = diffs[n-1]
    for k in range(n-2, -1, -1):
        s = s * (x - xs[k]) + diffs[k]

    return s
