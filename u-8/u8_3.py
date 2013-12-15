#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import exp, sin, cos, radians
import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt

def nth(arr, n):
    return arr[n]

def first(arr):
    return nth(arr, 0)
def second(arr):
    return nth(arr, 1)
def third(arr):
    return nth(arr, 2)
def forth(arr):
    return nth(arr, 3)

def runge_kutta(n, h, y0, F, t0=0):
    tn = float(t0)
    y = y0
    hh = h/2.0
    ys = []
    for i in range(1, n+1):
        K1 = F(tn, y)
        K2 = F(tn + hh, y + hh*K1)
        K3 = F(tn + hh, y + hh*K2)
        K4 = F(tn + h,  y + h*K3)
        y = y + (h/6.0)*(K1 + 2*K2 + 2*K3 + K4)
        ys.append(y)
        tn = tn + h

    return ys

def double_pendulum_f(t, y):
    theta1, theta2, omega1, omega2 = y
    g = 9.81
    dt = theta1 - theta2
    divisor = 2.0 - (cos(dt)**2)

    dividend1 = sum(
        [- omega1**2 * sin(dt) * cos(dt),
         + g * sin(theta2) * cos(dt),
         - omega2**2 * sin(dt),
         - 2.0 * g * sin(theta1)])

    dividend2 = sum(
        [(omega2**2 * sin(dt) * cos(dt)),
         + (2.0 * g * sin(theta1) * cos(dt)),
         + (2.0 * omega1**2 * sin(dt)),
         - (1.0 * g * sin(theta2))])

    return np.array([omega1, omega2, dividend1/divisor, dividend2/divisor])

## 8.3.3

def positions(ys):
    x1s = [cos(x) for x in map(first, ys)]
    x2s_raw = [cos(x) for x in map(second, ys)]
    x2s = [a + b for (a,b) in zip(x1s, x2s_raw)]

    y1s = [sin(x) for x in map(first, ys)]
    y2s_raw = [sin(x) for x in map(second, ys)]
    y2s = [a + b for (a,b) in zip(y1s, y2s_raw)]

    return x1s, y1s, x2s, y2s

y0 = np.array([radians(0.5), radians(0.5), 0.0, 0.0])

ys = runge_kutta(1000, 0.01, y0, double_pendulum_f)

# And now for the positions
_, _, x2s, y2s = positions(ys)
plt.plot(x2s, y2s)
plt.title('Position des 2. Pendels')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.show()

# 8.3.4

for (t1, t2) in [(0, 60), (160, 60)]:
    for (n, h) in [(500, 0.2), (1000, 0.1), (10000, 0.01)]:
        y0 = np.array([radians(t1), radians(t2), 0.0, 0.0])
        ys = runge_kutta(n, h, y0, double_pendulum_f)
        _, _, posx, posy = positions(ys)
        plt.plot(posx, posy, label=r'n = %d, h = %.3f' % (n, h))

    plt.title(r'Position des 2. Pendels bei $\theta_1 = %d^\circ,\ \theta2 = %d^\circ$' % (t1, t2))
    plt.legend()
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.show()
