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

y0 = np.array([radians(0.5), radians(0.5), 0.0, 0.0])

ys = runge_kutta(1000, 0.01, y0, double_pendulum_f)

plt.plot(range(1000), map(first, ys), label=r'$\theta_1$')
plt.plot(range(1000), map(second, ys), label=r'$\theta_2$')
plt.xlabel('Schritt')
plt.ylabel('Winkel')
plt.legend()
plt.show()
