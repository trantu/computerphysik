#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# 10.1.4

from math import exp
import numpy as np
import matplotlib.pyplot as plt

def u(x, t, U):
    if t == 0:
        return U/2.0 * exp(-abs(x)*U)

    raise NotImplementedError

def imp_euler(f, tN, dt, xs, dx):
    """Implicit euler for particle difussion"""

    xN = len(xs)
    mat = np.matrix(np.zeros((xN, tN)), dtype=float)
    r = dt / dx**2

    # Equivalent form for the implicit formula is
    # -r*u_{i-1}^n + (1+2r)*u_i^n - r*u_{i+1}^n = u_i^{n-1}

    B = np.array([(f(x, 0) for x in xs)])
    A = np.matrix(np.zeros((xN, xN)), dtype=float)

    # The first and last need to skip the first and last resp.
    A[0,0] = 1 + 2*r
    A[0,1] = -r
    A[-1,-2] = 1 + 2*r
    A[-1,-1] = -r

    # The rest is filled with [-r, 1+2r, -r]
    for i in range(1, xN-1):
        A[i,i-1] = -r
        A[i,i]   = 1 + 2*r
        A[i,i+1] = -r

    mat = np.array(np.zeros((tN, xN)), dtype=float)

    print 'B', B.shape, 'A', A.shape
    res = np.linalg.solve(A, B)
    print res[0]
    mat[0] = res
    for i in range(1, tN):
        B = res
        res = np.linalg.solve(A, B)
        mat[i] = res

    print mat
    return mat[:,-1].tolist()

tN = 6
dt = 0.1
dx = 0.1
xs = np.arange(-0.5, 0.51, dx)
res = imp_euler(lambda x, t: u(x, t, 10), tN, dt, xs, dx)
print res
[ys] = res
plt.plot(xs, ys)
plt.show()
