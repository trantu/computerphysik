#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

from numpy import matrix, zeros, eye
from math import sqrt
from functools import reduce

def alpha(a, pivot):
    return -1*(a/pivot)

def swap(m, a, b):
    temp = matrix(m[a])
    m[a] = m[b]
    m[b] = temp

# first, d is ignored
def gauss_alg(F, d):
    zeilen,spalten = F.shape
    # nicht moeglich zu eliminieren
    if zeilen != spalten or zeilen == 0 or spalten == 0:
        raise ValueError("ungleiche Anzahl zwischen Zeilen und Spalten")

    Ls = []

    # axis = 1 -> vertical
    #F = concatenate((F, d.T), axis=1)
    for s in range(spalten):
        for z in range(s,zeilen):
            pivot = F[z,z]

            #Zeilenvertauschen
            if pivot == 0:
                for ii in range(z+1,zeilen):
                    if F[ii,s] != 0:
                        swap(F, z, ii)
                        break # back to second loop

            pivot = F[z,z]
            for zz in range(z+1,zeilen):
                if F[zz,s] != 0:
                    alp = alpha(F[zz,s], pivot)
                    F[zz] = (alp * F[z]) + F[zz]
                    l = eye(zeilen, spalten)
                    l[zz] = (alp * l[z]) + l[zz]
                    Ls.append(l)

    # L3(L2(L1))
    L = reduce(lambda acc, x: x.dot(acc), Ls)

    return F, L

def backwards_substitution(R, y):
    xs = [0] * len(y)

    for i in reversed(range(len(y))):
        print "setting xs[%d] = %f" % (i, y[i,0])
        xs[i] = y[i,0]
        for j in reversed(range(i+1, len(y))):
            print "(%d,%d)" % (i, j)
            print "%f * %f" % (R[i,j], xs[j])
            xs[i] -= R[i,j] * xs[j]
        print "dividing xs[%d] (%f) / R[%d,%d] (%f)" % (i, xs[i], i, i, R[i,i])
        xs[i] = xs[i] / R[i,i]
    #xs.reverse()
    return xs

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

F, L = gauss_alg(matrix(A), None)

print "\nLR partition, L:\n", L

if (L.dot(A) == F).all():
    print "\nLR partition correct: LA = F ✔"
else:
    print "Oh no! LA != F"

R = L.dot(A)

b = matrix('-8; 15; 34', dtype=float)

# Transformed solution
Lb = L.dot(b)

xs = backwards_substitution(R, Lb)
print xs

print "control"

A = matrix('1 2 3; 6 -2 2; -3 1 -4', dtype=float)
F, L = gauss_alg(matrix(A, dtype=float), None)

print L
Lb = L.dot(matrix('12; -16; 2', dtype=float))

print "Lb\n", Lb

R = L.dot(A)
xs = backwards_substitution(R, Lb)
print xs
