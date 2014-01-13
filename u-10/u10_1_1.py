#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# 10.1.1
from math import exp,cos,sin
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

def u(x, t, U):
    if t == 0:
        return U/2.0 * exp(-abs(x)*U)

    raise NotImplementedError
def Bn(U,n):
    return quad(lambda x: U*np.exp(-abs(x)*U)*cos(2*np.pi*n*x), -0.5, 0.5)
def Cn(U,n):
    return quad(lambda x: -U**2*np.exp(-abs(x)*U)/(2*Bn(U,n)*np.pi*n), -0.5, 0.5)

def u10_1_1(f, ts, dt, xs, dx):
    U = 10
    NN = 100
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
            # since we're calculating u_i^{n+1}, we need to replace n -> n-1 in code
            for nn in range(NN):
                #mat[i,n] = mat[i,n-1] + r * (mat[i-1,n-1] - 2*mat[i,n-1] + mat[i+1,n-1])
                mat[i,n] = mat[i,n] + Bn(U,nn)*cos(2*np.pi*n*xs(i))*(Cn(U,nn)*sin(2*np.pi*n*ts(n))+cos(2*np.pi*n*ts(n)))
    return mat[:,-1].tolist()

# 10.1.2
t_lims = [0.0, 0.01, 0.03, 0.04, 0.05]

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
    ys = u10_1_1(lambda x, t: u(x, t, U), ts, dt, xs, dx)
    plt.plot(xs, ys, label='t=%g' % (t))

plt.title('Teilchendifussion , $dt=$%g, $dx=$%g, $r=$%g' % (dt, dx, r))
plt.ylabel('$u(x, t)$')
plt.xlabel('$x$')
plt.legend()
plt.show()
