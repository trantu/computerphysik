#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Carlos Martín, Tran Tu

import os, time
from math import sqrt, cos, pi
from itertools import islice
import matplotlib.pyplot as plt
import numpy as np

def take(seq, n):
    return list(islice(seq, 0, n))

class RandomGenerator:
    # Seed the random number generator
    def __init__(self, t=time.time()):
        # Start with the magic numbers
        pid = os.getpid()
        self.a = 7**5
        self.c = 0
        self.m = 2**31 - 1

        # And seed the RNG as per formula (2)
        num = abs((64979 * t * (pid - 83)) % 104729)
        num = num / 104729
        num = num * self.m
        self.num = int(num)

    # This makes this class an iterator itself
    def __iter__(self):
        return self

    # Sensible name and lets us make this an iterator
    def next(self):
        oldval = self.num

        # Calculate the new number, as per (1)
        self.num = (self.a * self.num + self.c) % self.m

        return oldval

# 12.1.2

def mean(xs):
    return sum(xs) / len(xs)

def sigma(xs):
    xmean = mean(xs)
    sq_diffs = [(xmean - x)**2 for x in xs]
    sq_mean = sum(sq_diffs) / len(xs)

    return sqrt(sq_mean)

m = 2**31 - 1
N = 10000

rng = RandomGenerator()
xs = [float(r) / m for r in take(rng, N)]


# squared differences from the mean

xmean = mean(xs)
xsig = sigma(xs)
print 'Mean', xmean
print 'σ   ', xsig

# Histogramm
plt.hist(xs, bins=np.arange(0, 1.01, 0.1))
plt.title(r'Histogramm der Zufallszahlen, Mean %f, $\sigma=$%f' % (xmean, xsig))
plt.xlabel('Werte')
plt.ylabel('Anzahl')
plt.show()

# 10.1.3

def C(xs, ys):
    if len(xs) != len(ys):
        raise ValueError

    xmean = mean(xs)
    xsig = sigma(xs)
    ymean = mean(ys)
    ysig = sigma(ys)
    n = len(xs)

    s = sum([(xs[i] - xmean) * (ys[i] - ymean) / (xsig * ysig) for i in range(n)])

    return s / n

ls = range(100, 10001, 100)
corrs = []
for l in ls:
    xs = take(rng, l)
    corrs.append(C(xs[:-1], xs[1:]))

plt.plot(ls, corrs, label='VL')

corrs = []
for l in ls:
    xs = take(rng, l)
    xcos = [cos((pi/180)*j + xj) for (xj, j) in zip(xs, range(1, len(xs)+1))]
    corrs.append(C(xcos[:-1], xcos[1:]))

plt.plot(ls, corrs, label='cos')

plt.title('Korrelation zwischen Zufallszahlen')
plt.legend()
plt.show()
