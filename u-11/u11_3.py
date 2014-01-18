#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Carlos Martín, Tran Tu

from math import sin, cos, pi, exp
from numpy.fft import fft, fftfreq, fftshift, ifft
import matplotlib.pyplot as plt
import random
import numpy as np

def s(t):
    return exp(-5 * (t - 1.5)**2) + 0.5 * exp(-2 * (t - 3)**4)

def n(t):
    return random.uniform(-0.5, 0.5)

def c(t):
    return s(t) + n(t)


random.seed(5)
dt = 0.02
ts = np.arange(0, 5, dt)
cs = map(c, ts)
ss = map(s, ts)

# 11.3.1
plt.title('Signal und Gesamtsignal')
plt.plot(ts, cs, label='Gesamtsignal')
plt.plot(ts, ss, label='Signal')
plt.legend()
plt.title('Signal')
plt.show()

# 11.3.2
freqs = fftfreq(len(ts), d=dt)
values = fft(cs)
ys = [abs(x)**2 for x in values]
plt.semilogy(fftshift(freqs), fftshift(ys))
plt.ylabel('Power')
plt.xlabel('Hz')
plt.title('Spektralanalyse')
plt.show()

# 11.3.3

# Wir wählen ω_0 5Hz. Da scheint es, keine große änderung zu geben,
# |N(x)|² = 90, aus dem Plot.

# The original says to use ω, but we don't that's a bit annoying to
# get. As |N(x)|² is a constant, we just need to match |C(w)|² with
# |C(w)|² which is obviously correct.
def phi(Ci, N = 90):
    C2 = abs(Ci)**2
    return (C2 - N) / C2

def filter_spectrum(C, i, w, w0, phi):
    if abs(w) < w0:
        return C[i] * phi(C[i], N=1)
    else:
        return 0

# Now that we have our filtering set up, let's take our fft values and
# filter them

C = fftshift(fft(cs))
freqs = fftshift(fftfreq(len(ts), d=dt))
w0 = 5

filtered = [filter_spectrum(C, i, freqs[i], w0, phi) for i in range(len(freqs))]
filtered_c = ifft(filtered)

plt.plot(ts, cs, label='Gesamtsignal')
plt.plot(ts, ss, label='Signal')
plt.plot(ts, [abs(x) for x in filtered_c], label='Filter')

plt.xlabel('Zeit')
plt.ylabel('Signal')
plt.legend()
plt.show()
