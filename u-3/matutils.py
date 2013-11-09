#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

from numpy import matrix

def alpha(a, pivot):
    return -1*(a/pivot)

def swap(m, a, b):
    temp = matrix(m[a])
    m[a] = m[b]
    m[b] = temp

# Perform Gauß triangulation while returning a generator of the steps it takes
def gauss_triangulation(A):
    # Make a copy so we can perform our calculations without disturting
    # the caller's data
    F = matrix(A)
    n, m = F.shape

    # nicht moeglich zu eliminieren
    if n != m or n == 0 or m == 0:
        raise ValueError("Check your matrix shape")

    for s in range(n):
        for z in range(s, m):
            pivot = F[z,z]

            #Zeilenvertauschen
            if pivot == 0:
                for ii in range(z+1, m):
                    if F[ii,s] != 0:
                        yield ('swap', z, ii)
                        swap(F, z, ii)
                        break # back to second loop

            pivot = F[z,z]
            for zz in range(z+1, m):
                if F[zz,s] != 0:
                    alp = alpha(F[zz,s], pivot)
                    yield ('div', alp, z, zz)
                    F[zz] = (alp * F[z]) + F[zz]

def backwards_substitution(R, y):
    """ Perform backwards subtitution on solution y and right-matrix R """
    xs = [0] * len(y)

    for i in reversed(range(len(y))):
        xs[i] = y[i,0]
        xs[i] -= sum([R[i,j] * xs[j] for j in range(i+1, len(y))])
        xs[i] /= R[i,i]

    return matrix(xs).transpose()
