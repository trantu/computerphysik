#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

import numpy as np
from numpy import matrix, zeros, eye, dot
from math import sqrt
from functools import reduce
from matutils import alpha, swap, gauss_triangulation, backwards_substitution
from u3_1 import gauss_alg

def flip(fn):
    return lambda a, b: fn(b, a)

def lr_partition(A):
    n, m = A.shape

    Ls = []

    for ins in gauss_triangulation(A):
        if ins[0] == 'div':
            _, alpha, a, b = ins
            l = matrix(eye(n,m), dtype=float)
            l[b] = (alpha * l[a]) + l[b]
            Ls.append(l)

    # L3(L2(L1))
    return reduce(flip(dot), Ls)

def cholesky(A):
    n, m = A.shape

    if n != m:
        raise AttributeError, "Matrix is not square"

    r = matrix(zeros((n, m)), dtype=float)

    for i in range(0, n):
        s = A[i,i] - sum([r[k,i]**2 for k in range(i)])

        if s < 0:
            raise AttributeError, "Matrix is not definite positive"
        
        r[i,i] = sqrt(s)
        for j in range(i+1, n):
            r[i,j] = (1.0/r[i,i]) * (A[i,j] - sum([r[k,i]*r[k,j] for k in range(i)]))

    return r

A = matrix([[64, -40, 16],
            [-40, 29, -4],
            [16,  -4, 62]], dtype=float)

CR = cholesky(A)

print "Cholesky R =\n", CR

if (CR.transpose().dot(CR) == A).all():
    print "\nCalculations correct: R^T * R = A ✔"
else:
    print "Oh no! R^T * R != A"

F = gauss_alg(matrix(A), None)
L = lr_partition(A)

print "\nLR partition, L:\n", L

if np.allclose(L.dot(A),F):
    print "\nLR partition correct: LA = F ✔"
else:
    print "Oh no! LA != F"

# Let's solve with b = [-8; 15; 34]

R = L.dot(A)
b = matrix('-8; 15; 34', dtype=float)

# Transformed solution
Lb = L.dot(b)

xs = backwards_substitution(R, Lb)
print "x = \n", xs

if np.allclose(dot(A, xs), b):
    print "\nLR partition worked: Ax = b ✔"
else:
    print "Oh no! Ax != b"

# print "\n*** Control: solution from the lecture ***"

# A = matrix('1 2 3; 6 -2 2; -3 1 -4', dtype=float)
# F, L = gauss_alg(matrix(A, dtype=float), None)

# print L
# Lb = L.dot(matrix('12; -16; 2', dtype=float))

# R = L.dot(A)
# xs = backwards_substitution(R, Lb)
# print "x = \n", xs
