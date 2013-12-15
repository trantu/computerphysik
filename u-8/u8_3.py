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
def fourth(arr):
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
    theta1s = map(first, ys)
    theta2s = map(second, ys)

    x1s = [sin(x) for x in theta1s]
    x2s = [x1 + sin(t1) for (x1, t1) in zip(x1s, theta1s)]

    y1s = [-cos(x) for x in theta1s]
    y2s = [y1 - cos(t2) for (y1, t2) in zip(y1s, theta2s)]

    return x1s, y1s, x2s, y2s


# And now for the positions
def auf_8_3_3():
    y0 = np.array([radians(0.5), radians(0.5), 0.0, 0.0])
    ys = runge_kutta(1000, 0.01, y0, double_pendulum_f)
    _, _, x2s, y2s = positions(ys)
    plt.plot(x2s, y2s)
    plt.title(r'Position 2. Pendel bei $\theta_1 = 0.5^\circ$, $\theta_2 = 0.5^\circ$')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.show()

# 8.3.4

def auf_8_3_4():
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

def auf_8_3_5():
    y0 = np.array([radians(0.5), radians(0.5), 0.0, 0.0])
    ys = runge_kutta(5000, 0.01, y0, double_pendulum_f)
    _, y1s, _, y2s = positions(ys)
    w1s = map(third, ys)
    w2s = map(fourth, ys)

    plt.suptitle('Winkelgeschwindigkeit gegen Position')
    plt.subplot(1, 2, 1)
    plt.plot(w1s, y1s, label='Pendel 1')
    plt.xlabel(r'$\omega_1$')
    plt.ylabel('$y_1$')

    plt.subplot(1, 2, 2)
    plt.plot(w2s, y2s, label='Pendel 2')
    plt.xlabel(r'$\omega_2$')
    plt.ylabel('$y_2$')

    plt.show()

# 8.3.6
def auf_8_3_6():
    y0 = np.array([radians(160), radians(60), 0.0, 0.0])
    ys = runge_kutta(1000000, 0.01, y0, double_pendulum_f)
    #ys = runge_kutta(100000, 0.01, y0, double_pendulum_f)

    w1s = map(third, ys)
    w2s = map(fourth, ys)

    H, xedges, yedges = np.histogram2d(w1s, w2s, bins=[256, 256])
    plt.imshow(H, extent=[yedges[0], yedges[-1], xedges[-1], xedges[0]])
    plt.xlabel(r'$\omega_1$')
    plt.ylabel(r'$\omega_2$')
    plt.title('Histogramm der Wilkelgeschwindigkeiten')
    plt.show()

#auf_8_3_3()
#auf_8_3_4()
#auf_8_3_5()
#auf_8_3_6()

