#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Carlos Martín, Tran Tu

from math import sin, cos, pi
from numpy.fft import fft, fftfreq, fftshift
import matplotlib.pyplot as plt
import numpy as np

def signal(x):
    return sin(x)  + 2 * sin(0.6*x) + 0.5 * sin(4*x)

def blackman(i, N):
    piN = (pi * i) / N
    return 0.42 - 0.5 * cos(2*piN) + 0.08 * cos(4*piN)

def trapez(i, N):
    if i < N/10.0:
        return i / (N/10.0)
    elif i > (N - N/10.0):
        adj_i = i - (N - N/10.0) # adjusted i, starting at top 10%
        return 1 - adj_i / (N/10.0)
    else:
        return 1.0

window_fns = [
    (lambda i, N: 1, 'Straight'),
    (blackman, 'Blackman'),
    (trapez, 'Trapez'),
    ]

for N in [2000, 10000]:
    plt.figure()
    plt.title('Spektralanalyse N = %d' % N)
    freqs = fftfreq(N, d=0.05)
    xs = [0.05 * i for i in range(N)]
    a = map(signal, xs)

    for (w, lbl) in window_fns:
        in_values = [a[i] * w(i, N) for i in range(N)]
        values = [abs(x)**2 for x in fft(in_values)]
        plt.semilogy(fftshift(freqs), fftshift(values), label=lbl)

        plt.legend()

plt.show()
