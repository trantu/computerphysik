#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from math import exp
import numpy as np
import matplotlib.pyplot as plt

def u(x, t, U):
    if t == 0:
        return U/2.0 * exp(-abs(x)*U)

    raise NotImplementedError

# 10.1.7
def crank_nicolson(f, ts, dt, xs, dx):
    """Crank-Nicolson for particle difussion"""

    xN = len(xs)
    tN = len(ts)
    mat = np.matrix(np.zeros((xN, tN)), dtype=float)
    r = dt / dx**2

    # Equvalent form of the Crank-Nicolson multiplied with dt
    # -(r/2)u_{i-1}^n + (1+r)u_i^n -(r/2)u_{i+1}^n =
    # = (r/2)u_{i-1}^{n-1} + (1-r)u_{i}^{n-1} + (r/2)u_{i+1}^{n-1}

    B = np.array([f(x, 0) for x in xs])
    A = np.matrix(np.zeros((xN, xN)), dtype=float)

    # The first and last need to skip the first and last resp.
    A[0,0] = 1 + r
    A[0,1] = -r/2.0
    A[-1,-2] = -r/2.0
    A[-1,-1] = 1 + r

    # The rest is filled with [-r/2, 1+r, -r/2]
    for i in range(1, xN-1):
        A[i,i-1] = -r/2.0
        A[i,i]   = 1 + r
        A[i,i+1] = -r/2.0

    mat = np.array(np.zeros((tN, xN)), dtype=float)

    res = np.linalg.solve(A, B)
    mat[0] = res
    # Now we loop and create the new B based off of the results of the
    # last solution.
    for n in range(1, tN):
        B = np.array(np.zeros((xN,)), dtype=float)

        # The first and last are different, as usual
        B[0] = (1-r)*res[0] + (r/2.0)*res[1]
        B[-1] = (r/2.0)*res[-2] + (r-1)*res[-1]

        # The rest is filled with r/2 + (1-r) + r/2
        for i in range(1, xN-1):
            B[i] = (r/2.0)*res[i-1] + (1-r)*res[i] + (r/2.0)*res[i+1]

        res = np.linalg.solve(A, B)
        mat[n] = res

    return mat[-1,:].tolist()

U = 10
tN = 6
dx = 0.01
t_lims = [0.0, 0.01, 0.03, 0.04, 0.05]

for dt in [0.01, 1e-5]:
    plt.figure()
    for t in t_lims:
        xs = np.arange(-0.5, 0.5+dx, dx)
        ts = np.arange(0, t+dt, dt)
        ys = crank_nicolson(lambda x, t: u(x, t, U), ts, dt, xs, dx)
        plt.plot(xs, ys, label='t=%g' % (t))

        r = dt / dx**2
        plt.title('Teilchendifussion (Crank-Nicolson) $r=$%g' % (r))
        plt.xlabel('$x$')
        plt.ylabel('$u(x, t)$')
        plt.legend()

plt.show()

# Ähnlich zu den anderen Verfahren, konvergiert der Verfahren abhängig
# von r. Im ersten Graph ist r=100 und der Verfahren konvergiert
# nicht. Im zwiten Graph ist r=0.1 und der Verfahren konvergiert.
