#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

from numpy import matrix, zeros
from math import sqrt

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

R = cholesky(A)

print "R =\n", R

if (R.transpose().dot(R) == A).all():
    print "\nCalculations correct: R^T * R = A ✔"
else:
    print "Oh no! R^T * R != A"
