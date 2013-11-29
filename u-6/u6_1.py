#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tu Tran

import numpy as np
import matplotlib.pyplot as plt
from math import cos

def ansatz1():
    def f0(x):
        return 1

    def f1(x):
        return (x - 1970.0) / 100.0

    def f2(x):
        return ((x - 1970.0) / 100.0)**2

    return [f0, f1, f2]

def ansatz2():
    def f0(x):
        return cos((x - 1970.0) / 100.0)

    def f1(x):
        return (x - 1970.0) / 100.0

    def f2(x):
        return 0

    def f3(x):
        return ((x - 1970.0) / 100.0)**3

    # return [f0, f1, f2, f3]
    return [f0, f1, f3]

def curve_fitting(xs, ys, fs):
    """Fit the curve given by data and fs. data is an array with x-y pairs
    (e.g. from numpy.recfromtxt()) and fs is a list of f_i
    """

    A = []
    ys = np.matrix(ys).transpose()
    for x in xs:
        row = [f(x) for f in fs]
        A.append(row)
    A = np.matrix(A)

    AT = A.transpose()
    ATA = np.dot(AT, A)
    ATy = np.dot(AT, ys)

    a = np.linalg.solve(ATA, ATy)

    assert len(fs) == len(a)

    def fit(x):
        return sum([a[i,0] * fs[i](x) for i in range(len(fs))])

    return a, fit

def show_errors(xs, ys, f):
    errors = [(xs[i], abs(f(xs[i]) - ys[i])) for i in range(len(xs))]
    maxerr = max(errors, key=lambda (_, e): e)
    errorsq = sum([x*x for (_, x) in errors])

    print u'Fehlerquadrate: ', errorsq
    print u'Maximaler fehler ist bei %.1f -> %f ' % maxerr
    

data = np.recfromtxt('data.txt')
xs = data[:,0]
ys = data[:,1]
a, fit1 = curve_fitting(xs, ys, ansatz1())
print "### Funktion i)"
print "Optimale Parameter"
print(a)
print ""
print u'Dönerpreis in 2020:', fit1(2020)
show_errors(xs, ys, fit1)

a, fit2 = curve_fitting(xs, ys, ansatz2())
print "### Funktion ii)"
print "Optimale Parameter"
print(a)
print ""
print u'Dönerpreis in 2020:', fit2(2020)
show_errors(xs, ys, fit2)

plt.plot(xs, ys, label=ur'Original')

points = np.arange(1970, 2010, 0.1)
plt.plot(points, map(fit1, points), label=ur'fit1')
plt.plot(points, map(fit2, points), label=ur'fit2')
plt.legend(loc='upper left')
plt.show()
