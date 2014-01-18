#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Carlos Martín, Tran Tu

from math import sin, cos, pi, exp
from numpy.fft import fft, fftfreq
import matplotlib.pyplot as plt
import random
import numpy as np

def s(t):
    return exp(-5 * (t - 1.5)**2) + 0.5 * exp(-2 * (t - 3)**4)

def n(t):
    return random.uniform(-0.5, 0.5)

def c(t):
    return s(t) + n(t)


# 11.3.1
random.seed(5)
dt = 0.02
ts = np.arange(0, 5, dt)

plt.plot(ts, map(s, ts), label='Signal')
plt.plot(ts, map(c, ts), label='Gesamtsignal')
plt.show()

# 11.3.2
freqs = fftfreq(len(ts), d=dt)
values = fft(map(c, ts))
ys = [abs(x)**2 for x in values]
plt.semilogy(freqs, ys)
plt.show()
