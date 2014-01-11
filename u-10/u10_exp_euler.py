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

def exp_euler(f, tN, dt, xs, dx):
    """Explicit Euler for particle diffussion
    f: function to use
    tN: how many times to use
    dt: difference in time
    xs: list of x values
    dx: difference in space
    """

    xN = len(xs)
    mat = np.matrix(np.zeros((xN, tN)), dtype=float)
    r = dt / dx**2

    # Initial values with u(x, t=0)
    for i in range(len(xs)):
        mat[i,0] = f(xs[i], 0)
    # set u(∓1/2,t) = 0
    for n in range(tN):
        mat[0,n] = 0

    for n in range(1, tN):
        for i in range(xN-1):
            # since we're calculating u_i^{n+1}, we need to replace n -> n-1 in code
            mat[i,n] = mat[i,n-1] + r * (mat[i,n-1] - 2*mat[i,n-1] + mat[i+1,n-1])

    return mat[:,-1].tolist()

# 10.1.4
ts = [0, 0.01, 0.03, 0.04, 0.05]

U = 10
tN = 6
dt = 0.01
dx = 0.1

#for dx in [0.1]:
for dx in [0.11, 0.15, 0.2, 0.3, 0.4]:
    xs = np.arange(-0.5, 0.5, dx)
    ys = exp_euler(lambda x, t: u(x, t, U), tN, dt, xs, dx)
    r = dt / dx**2
    plt.plot(xs, ys, label='dx = %f, r = %f' % (dx, r))

plt.ylabel('$u(x, t_{\text{end}})$')
plt.xlabel('$x$')
plt.legend()
plt.show()
