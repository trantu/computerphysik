#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Carlos Martín, Tran Tu

import os, time
from math import sqrt
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

m = 2**31 - 1
N = 10000

rng = RandomGenerator()
xs = [float(r) / m for r in take(rng, N)]
mean = sum(xs) / N

# squared differences from the mean
sq_diffs = [(mean - x)**2 for x in xs]
sq_mean = sum(sq_diffs) / N
stddev = sqrt(sq_mean)

print 'Mean', mean
print 'σ   ', stddev

# Histogramm
plt.hist(xs, bins=np.arange(0, 1.01, 0.1))
plt.title(r'Histogram der Zufallszahlen, Mean %f, $\sigma=$%f' % (mean, stddev))
plt.xlabel('Werte')
plt.ylabel('Anzahl')
plt.show()
