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

# 10.1.4
def imp_euler(f, ts, dt, xs, dx):
    """Implicit euler for particle difussion"""

    xN = len(xs)
    tN = len(ts)
    mat = np.matrix(np.zeros((xN, tN)), dtype=float)
    r = dt / dx**2

    # Equivalent form for the implicit formula is
    # -r*u_{i-1}^n + (1+2r)*u_i^n - r*u_{i+1}^n = u_i^{n-1}

    B = np.array([f(x, 0) for x in xs])
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

    res = np.linalg.solve(A, B)
    mat[0] = res
    for i in range(1, tN):
        B = res
        res = np.linalg.solve(A, B)
        mat[i] = res

    return mat[-1,:].tolist()

U = 10
tN = 6
dt = 0.01
dx = 0.1
r = dt / dx**2
t_lims = [0.0, 0.01, 0.03, 0.04, 0.05]

for t in t_lims:
    xs = np.arange(-0.5, 0.5+dx, dx)
    ts = np.arange(0, t+dt, dt)
    ys = imp_euler(lambda x, t: u(x, t, U), ts, dt, xs, dx)
    plt.plot(xs, ys, label='t=%g' % (t))

plt.title('Teilchendifussion (imp) $r=$%g' % (r))
plt.xlabel('$x$')
plt.ylabel('$u(x, t)$')
plt.legend()
plt.show()

# 10.1.5: Im gegesatz zum expliziten Euler-Verfahren, muss r hier r >1
# bleiben. Ändert man dx und dt (solange r > 1)und plottet man die
# Ergebnisse merkt man kein großes Unterschied bei der Plots. Daraus
# kann man die Schlussfolgerung ziehen, dass diese Verfahren ziemlich
# stabil ist.

# 10.1.6: Jeder Schritt beim Impliziten Euler-Verfahren ist teurer als
# beim expliziten (LGS lösen gegen Multiplikationen und Additionen in
# einer Schleife). Man würde denken, dass das explizite Verfahren
# billiger sein müsste.
#
# Es gibt aber ein wichtiger Unterschied, und zwar was es 'r'
# angeht. Beim expliziten Verfahren muss r<1 sein, was bedeutet, dass
# dt ein kleiner Wert ist, was wiederum bedeutet, dass wir viele
# Durchläufe der inneren Schleife machen müssen, bis wir zum
# gewünschten 't' ankommen.
#
# Beim impliziten Verfahren ist zwar jeder Schritt teurer, aber da r>1
# sein muss, ist dt deutlich größer als beim expliziten Verfahren. In
# Größeres dt bedeutet viele weniger Durchläufe. Bei gleichem dt wäre
# dieses langsamer, aber bei dieser Aufgabe ist das implizite
# Verfahren deutlich schneller, da wir dt = 1e-2 statt dt = 1e-5
# benutzen können.
