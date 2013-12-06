#!/usr/bin/evn python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

import math
from math import sin, cos
from numpy import linspace
import matplotlib.pyplot as plt

ur"""Zweite ableitung

\begin{align*}
f''(x) &= \frac{1}{2h} ( Df(x+h, h) - Df(x-h, h )\\
&= \frac{1}{2h} \left[ \frac{f(x+h+h) - f(x+h-h}{2h} + O(h^2) -
\frac{f(x-h+h) - f(x-h-h)}{2h} + O(h^2) \right] + O(h^2)\\
&= \frac{1}{2h}\left[ \frac{1}{2h}\left[ f(x+2h) - f(x) - f(x) + f(x-2h)\right] + O(h^2) + O(h^2) \right] + O(h^2)\\
&= \frac{f(x+2h) - 2f(x) + f(x-2h)}{4h^2} + O(h^2)
\end{align*}

"""

# did someone say lisp?
def second((_, x)):
    return x

def f(x):
    return 1.0 / (2 + cos(x))

def f_1(x):
    return sin(x) / (cos(x) + 2)**2

def Df(x, h):
    return (f(x+h) - f(x-h)) / (2*h)

def D1f(x, h):
    return (f(x+h) - f(x)) / h

xs = linspace(0, 2*math.pi, 50)
hs = [(2**(-i), i) for i in range(6)]

plt.figure()
# This one doesn't change
plt.plot(xs, map(f_1, xs), label=r"$f'(x)$")
for (h, i) in hs:
    plt.plot(xs, [D1f(x, h) for x in xs], label=u'i = %d' % i)

plt.title('Ableitung durch $D_1f$')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='lower left')

plt.figure()
plt.plot(xs, map(f_1, xs), label=r"$f'(x)$")
for (h, i) in hs:
    plt.plot(xs, [Df(x, h) for x in xs], label=u'i = %d' % i)

plt.title('Ableitung durch $Df$')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='lower left')

## h-Extrapolation

def richardson(x, h, n, fn=Df):
    D0 = [fn(x, h/2.**i) for i in range(n+1)]
    
    def D(i, k):
        if k == 0:
            return D0[i]

        tmp = D(i+1, k-1)
        return tmp + (tmp - D(i, k-1)) / (2.**k - 1)

    return D(0, n)

hs = [(2**(-k), k) for k in range(11)]

plt.figure()
err = []
for (h, k) in hs:
    err.append(sum([abs(Df(x, h) - f_1(x)) for x in xs]))
plt.semilogy(map(second, hs), err, label='$Df$')

err = []
for (h, k) in hs:
    err.append(sum([abs(D1f(x, h) - f_1(x)) for x in xs]))
plt.semilogy(map(second, hs), err, label='$D_1f$')

err = []
for (h, k) in hs:
    err.append(sum([abs(richardson(x, h, 2, fn=Df) - f_1(x)) for x in xs]))
plt.semilogy(map(second, hs), err, label='Richardson $Df$')

err = []
for (h, k) in hs:
    err.append(sum([abs(richardson(x, h, 1, fn=D1f) - f_1(x)) for x in xs]))
plt.semilogy(map(second, hs), err, label='Richardson $D_1f$')

plt.ylabel('Fehler')
plt.xlabel('k')
plt.legend(loc='lower left')

# Show the plots
plt.show()
