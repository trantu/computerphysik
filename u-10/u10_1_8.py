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
# Aufgabe 10.1.8
def exp_euler(f, ts, xs, dx):
    """Explicit Euler for particle diffussion
    f: function to use
    tN: how many times to use
    dt: difference in time
    xs: list of x values
    dx: difference in space
    """

    xN = len(xs)
    tN = len(ts)
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
            #TODO 
            # since we're calculating u_i^{n+1}, we need to replace n -> n-1 in code
            mat[i,n] = mat[i,n-1] + r * (mat[i-1,n-1] - 2*mat[i,n-1] + mat[i+1,n-1])-dt*5000*mat[i,n-1]*np.exp(-1*abs(xs[i]-0.5)*10)
    xK,xP = mittl_Konz_Pos(mat, tN, xN)
    return mat[:,-1].tolist(),xK,xP

def mittl_Konz_Pos(mat,tN,xN):
    xK = np.zeros(tN)
    xP = np.zeros(tN)
    for n in range(tN):
        for i in range(xN):
            xK[n] = xK[n]+mat[i,n]
            xP[n] = xP[n] + xs[i] * mat[i,n]
        xK[n] = xK[n]/xN
        xP[n] = xP[n]/xN
    return xK,xP
# 10.1.8
t_lims = np.arange(0, 0.06, 0.01)

U = 10
dt = 0.00001
dx = 0.01
r = dt / dx**2

# We should be able to do this calculation only once and look into
# different offsets in a returned matrix, but this works, and isn't
# overly expensive

for t in t_lims:
    xs = np.arange(-0.5, 0.5+dx, dx)
    ts = np.arange(0, t+dt, dt)
    ys,_,_ = exp_euler(lambda x, t: u(x, t, U), ts, xs, dx)
    plt.plot(xs, ys, label='t=%g' % (t))

xs = np.arange(-0.5, 0.5+dx, dx)
t_lims = np.arange(0, 0.2, 0.00001)    
ys,xK,xP = exp_euler(lambda x, t: u(x, t, U), t_lims, xs, dx)
plt.plot(t_lims, xK)
plt.show()
plt.plot(t_lims, xP)
plt.show()
plt.title('Teilchendiffussion, $dt=$%g, $dx=$%g, $r=$%g' % (dt, dx, r))
plt.ylabel('$u(x, t)$')
plt.xlabel('$x$')
plt.legend()
plt.show()

