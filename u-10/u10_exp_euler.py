#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# 10.1.4

from math import exp
import numpy as np
import matplotlib.pyplot as plt

def u(x, t, U):
    if t == 0:
        return U/2.0 * exp(-abs(x)*U)

def exp_euler(f, tN, dt, xN, dx):
    """Explicit Euler for particle diffussion
    f: function to use
    dt: difference in time
    dx: difference in space
    """

    mat = np.matrix(np.zeros((tN+1, xN)), dtype=float)
    r = dt / dx**2

    # Initial values with u(x, t=0)
    for n in range(xN):
        mat[0,n] = f(n * dx, 0)

    print mat

    for i in range(1, tN):
        # copy over the results from the last
        mat[i,0] = mat[i-1,1]
        mat[i,xN-1] = mat[i-1, xN-1]
        for n in range(1, xN-1):
            # since we're calculating u_i^{n+1}, we need to replace n -> n-1 in code
            mat[i,n] = mat[i,n-1] + r * (mat[i,n-1] - 2*mat[i,n-1] + mat[i+1,n-1])

    return mat[-2,:]

# 10.1.4
ts = [0, 0.01, 0.03, 0.04, 0.05]

U = 10
tN = 6
dt = 0.01
xN = 15
dx = 1

[ys] = exp_euler(lambda x, t: u(x, t, U), tN, dt, xN, dx).tolist()
print ys
r = dt / dx**2
plt.plot(ys, label='r = %f' % r)
plt.legend()
plt.show()
