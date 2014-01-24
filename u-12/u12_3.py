#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Carlos Martín, Tran Tu

import os, time
from itertools import islice

def take(seq, n):
    return list(islice(seq, 0, n))

class RandomGenerator(object):
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

class UniformRNG(RandomGenerator):
    def next(self):
        n = super(UniformRNG, self).next()
        return n / float(self.m)

import itertools
import numpy as np
import matplotlib.pyplot as plt
from math import exp, sqrt, pi

N = 2
d = 2

def f(xs):
    return exp(-sum(xs**2))

def mittelpunkt(f, a, b, N, d):
    _a = a
    _b = b
    is_ = np.array(list(itertools.product(range(1, N+1), repeat=d)), dtype=float)
    a = np.array([a]*d, dtype=float)
    b = np.array([b]*d, dtype=float)

    tmp = sum([f(a + (b-a)*(i-0.5)/N) for i in is_])
    tmp *= (_b - _a) / float(N**d)

    return tmp

def hit_miss_monte(f, fmin, fmax, N, d, a, b):
    _a = a
    _b = b

    a = np.array([a]*d, dtype=float)
    b = np.array([b]*d, dtype=float)

    N_plus = 0.0
    rng = UniformRNG()

    for _ in range(N**d):
        ri = float(rng.next())
        si = float(rng.next())
        if f(a+(b-a)*ri) - fmin > (fmax-fmin)*si:
            N_plus = N_plus + 1

    return (_b - _a) * ((N_plus/float(N**d))*(fmax - fmin) + fmin)

def stuetz_monteC_integr(f, N, d, a, b):
    vorfaktor = (b - a)/float(N**d)
    a = np.array([a]*d, dtype=float)
    b = np.array([b]*d, dtype=float)

    summe = 0.0
    rng = UniformRNG()
    for _ in range(N**d):
        r = float(rng.next())
        summe = summe + f(a + (b-a) * r)

    return vorfaktor * summe

import scipy.special

def exact(d):
    return ((sqrt(pi)/2) * scipy.special.erf(1))**d

def calc_errors(N, ds):
    mitt_errs = []
    hm_errs = []
    stu_errs = []
    for d in ds:
        mitt = mittelpunkt(f, a=0, b=1, N=N, d=d)
        hm = hit_miss_monte(f, fmin=0, fmax=1, N=N, d=d, a=0, b=1)
        stu = stuetz_monteC_integr(f, N=N, d=d, a=0, b=1)
        ex = exact(d)

        err = abs((mitt - ex)/ex)
        mitt_errs.append((N, d, err))

        err = abs((hm - ex)/ex)
        hm_errs.append((N, d, err))

        err = abs((stu - ex)/ex)
        stu_errs.append((N, d, err))

    return mitt_errs, hm_errs, stu_errs

def plot_errors(mitt_errs, hm_errs, stu_errs, N):
    plt.figure(0) # for the semilog plot

    xs, ys = zip(*[(d, err) for N, d, err in mitt_errs])
    plt.semilogy(xs, ys, label='Mittelpunkt, $N=%d$' % N)

    xs, ys = zip(*[(d, err) for N, d, err in hm_errs])
    plt.semilogy(xs, ys, label='Monte Carlo, $N=%d$' % N)

    xs, ys = zip(*[(d, err) for N, d, err in stu_errs])
    plt.semilogy(xs, ys, label=u'Stützstellen Monte Carlo, $N=%d$' % N)

    plt.figure(1) # for the loglog plot

    xs, ys = zip(*[(N**d, err) for N, d, err in mitt_errs])
    plt.loglog(xs, ys, label='Mittelpunkt, $N=%d$' % N)

    xs, ys = zip(*[(N**d, err) for N, d, err in hm_errs])
    plt.loglog(xs, ys, label='Monte Carlo, $N=%d$' % N)

    xs, ys = zip(*[(N**d, err) for N, d, err in stu_errs])
    plt.loglog(xs, ys, label=u'Stützstellen Monte Carlo, $N=%d$' % N)


plot_errors(*calc_errors(2, range(1, 18)), N=2)
plot_errors(*calc_errors(10, range(1, 5)), N=10)

plt.figure(0)
plt.title('Relative Fehler bei der Integration, semilog')
plt.xlabel(r'$d$')
plt.ylabel('Relativer Fehler')
plt.legend(loc='upper left')

plt.figure(1)
plt.title('Relative Fehler bei der Integration, loglog')
plt.xlabel(r'$N^d$')
plt.ylabel('Relativer Fehler')
plt.legend(loc='upper left')

plt.show()
