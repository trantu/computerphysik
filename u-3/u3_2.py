#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

from numpy import matrix, zeros, eye
from math import sqrt

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

    L = []

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
                    #print F[z,:]
                    #print F[z+1,:]
                    alp = alpha(F[zz,s], pivot)
                    F[zz] = (alp * F[z]) + F[zz]
                    l = eye(zeilen, spalten)
                    l[zz] = (alp * l[z]) + l[zz]
                    L.append(l)

    LL = reduce(lambda acc, x: x.dot(acc), L)

    return F, LL

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

F, L = gauss_alg(matrix(A), None)

print F, L

# print "control"

# F, L = gauss_alg(matrix('1 2 3; 6 -2 2; -3 1 -4', dtype=float), None)

# print L
