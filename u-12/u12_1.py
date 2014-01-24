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
        # Calculate the new number, as per (1)
        self.num = (self.a * self.num + self.c) % self.m

        return self.num

# 12.1.2

def mean(xs):
    return sum(xs) / len(xs)

def sigma(xs):
    xmean = mean(xs)
    sq_diffs = [(xmean - x)**2 for x in xs]
    sq_mean = sum(sq_diffs) / len(xs)

    return sq_mean
    #return sqrt(sq_mean)

def u12_1_2():
    m = 2**31 - 1
    N = 10000

    rng = RandomGenerator()
    xs = [float(r) / m for r in take(rng, N)]


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

    return s / float(n)

ls = range(100, 10001, 100)

def u12_1_3():
    rng = RandomGenerator()
    corrs = []
    for l in ls:
        xs = take(rng, l)
        corrs.append(C(xs[:-1], xs[1:]))

    plt.plot(ls, corrs, label='VL')

    corrs = []
    for l in ls:
        xs = take(rng, l)
        xcos = [cos((pi/180)*j + xj) for (xj, j) in zip(xs, range(len(xs)))]
        corrs.append(C(xcos[:-1], xcos[1:]))

    plt.plot(ls, corrs, label='cos')

    plt.title('Korrelation zwischen Zufallszahlen')
    plt.legend()
    plt.show()

class DeterminateRandomGenerator(RandomGenerator):
    def __init__(self, a, c, m, x0=757685724):
        self.a = a
        self.c = c
        self.m = m
        self.num = x0

from mpl_toolkits.mplot3d import Axes3D

def u12_1_4():
    k = 1
    N = 10000 + 2
    magics = [(1,1,1024), (61,0,1024), (61,1,1024), (101,2,16384), (65539,0,2**31)]

    for args in magics:
        # We generate exta numbers so we can have all three lists
        # with 10000 numbers, the extra at the end are thrown away
        rng = DeterminateRandomGenerator(*args)
        nums = take(rng, N)
        xs = nums[:-2]
        ys = nums[1:-1]
        zs = nums[2:]

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xs, ys, zs, '.', s=1)
        plt.title('$a=$ %d, $c=$ %d, $m=$ %d' % args)
        plt.xlabel('x')
        plt.ylabel('y')

    plt.show()

if __name__ == '__main__':
    #u12_1_2()
    #u12_1_3()
    u12_1_4()
