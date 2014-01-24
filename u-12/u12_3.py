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
from math import exp

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
    tmp *= (_b - _a) / float(N)

    return tmp

if __name__ == '__main__':
    print mittelpunkt(f, a=0, b=1, N=10, d=1)
