#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import cos, pi, sqrt, log10, ceil

from itertools import cycle, islice, repeat, izip
from operator import sub

from numpy import linspace, matrix, zeros
import matplotlib.pyplot as plt

from scipy.special.orthogonal import p_roots

import numpy as np
from numpy.core.numeric import dtype

def take(seq, n):
    return list(islice(seq, 0, n))

def avsub(seq):
    return [abs(sub(*x)) for x in seq]

def f(x):
    return 1.0 / (2.0 + cos(x))

def rechteckregel(xs, a, b, f):
    n = len(xs)
    h = (b - a)/float(n)
    ys = [f(xs[i] + h/2.0) for i in range(n)]
    return h*sum(ys)

def trapezregel(xs, a, b, f):
    n = len(xs)
    h = (b - a)/float(n - 1)
    s = 0
    for i in range(1, n):
        s += f(xs[i])

    return h * ((f(a) - f(b)) / 2.0 + s)

# def simpsonregel(xs, a, b, f):
#     n = len(xs)
#     h = (b - a)/float(n)

#     return (h/3.0) * (f(b) + f(a) + 4*f(a+h))

def simpsonregel(xs, a, b, f):
    n = len(xs)
    h = (b - a)/float(n-1) # need the number of intervals

    coefficients = [1] + take(cycle([4,2]), n-2) + [1]
    lst = zip(coefficients, xs)
    s = map(lambda (c, x): c * f(x), lst)
    return (h/3.0) * sum(s)
# Aufgabe 7.2.4    
def gauss_quadratur(xs,a,b,f):
    
    [x,w] = p_roots(len(xs))
    summe = 0.0
    #print len(x)
    for i in range (0,len(x)):
        xu = (b-a)/2 * x[i] + (a+b)/2
        summe = summe + (w[i] * f(xu))
    return (b-a)/2 * summe
    """
    Bemerkung: 
    
    Die Gauß-Quadratur, wie das Verfahren durchgefuehrt wird, ist recht aufwendig, hat aber eine hohe Genauigkeit.
    Die Nullstellen des Legendre Polynoms werden numerisch bestimmt, dadurch sinkt die Effizienz. 
    Ein unlineares Gleichungssystem muss gelöst werden, damit optimale xi und ai bestimmt werden. 
    Falls f(x) an beliebigen Stellen verfuegbar ist, sind die Gauss-Quadraturformeln im allgemeinen 
    viel effizienter als z.B. das Romberbergg-Verfahren. Fuer die Fehlerschaetzung muss man den Aufwand vergroeßern.
    """
# Fake rule to show the right thing
def correct(_xs, _a, b, _f):
    if b == pi/2.0:
        return pi/(3*sqrt(3))
    if b == pi:
        return pi/sqrt(3)

    raise ValueError('Invalid b')

## Let's calculate some integrals
a, b, n = 0, pi/2.0, 21
xs = linspace(a, b, n)
funs = [
    ('Genau', correct),
    ('Rechteck', rechteckregel),
    ('Trapez', trapezregel),
    ('Simpson', simpsonregel),
    ]

print '### x0 = pi/2'
for (lbl, fn) in funs:
    print '%s: %f' % (lbl, fn(xs, a, b, f))
#print 'Gauss-Quadratur', gauss_quadratur(a, b, f, n)

a, b, n = 0, pi, 21
xs = linspace(a, b, n)
print '### x0 = pi'
for (lbl, fn) in funs:
    print '%s: %f' % (lbl, fn(xs, a, b, f))

def romberg(a, b, n):
    I = matrix(zeros((n, n)), dtype=float)
    # Fill in the initial column
    for i in range(n):
        h = (b / 4.0) * 2**(-i)
        xs = np.arange(a, b, h)
        I[i,0] = trapezregel(xs, a, b, f)

    # And now we iterate, column-first
    for k in range(1, n):
        for i in range(k, n):
            I[i,k] = I[i,k-1] + (I[i,k-1] - I[i-1,k-1]) / (4**k - 1)

    return I

# Let's show how Ii,1 and Ii,i differ from the real value

# have a list of the real solution at the ready
truth = repeat(pi/(3.0 * sqrt(3.0)))

a, b, n = 0, pi/2.0, 20
I = romberg(a, b, n)

#Ii,i
diffs = zip([I[i,i] for i in range(n)], truth)
plt.semilogy(range(n), avsub(diffs), label='$I_{i,i}$')

# I1,i
diffs = zip([I[1,i] for i in range(1, n)], truth)
plt.semilogy(range(1, n), avsub(diffs), label="$I_{1,i}$")

# Simpson
diffs = zip([simpsonregel(linspace(a, b, i+1), a, b, f) for i in range(2, n, 2)], truth)
plt.semilogy(range(2, n, 2), avsub(diffs), label="Simpson")

#plt.plot(range(5), [I[i,i] for i in range(5)])
#plt.plot(range(1, 5), [I[i,1] for i in range(1, 5)])
plt.legend(loc='lower left')
plt.show()

#print 'Gauss-Quadratur', gauss_quadratur(a, b, f, n)
#Aufgabe 7.2.5

funs = [
    ('Trapez', trapezregel),
    ('Simpson', simpsonregel),
    (u'Gauß', gauss_quadratur),
    ]

def a7_2_5():
    N = range(2,41,2)

    I0_1 = round(pi/(3*sqrt(3)),14)
    I0_2 = round(pi/sqrt(3),13)
    x0_1= pi/2.0
    x0_2 = pi
    for (lbl, fn) in funs:
        E = []
        for i in N:
            xs = linspace(0, x0_1, i+1)
            fehler=abs(fn(xs, 0, x0_1, f) - I0_1)
            E.append(fehler)
        plt.semilogy(N, E, label=lbl)

    I = romberg(0, x0_1, 4)
    plt.semilogy([4, 8, 16, 32], [abs(I[i,i] - I0_1) for i in range(4)], label='Romberg')
    plt.legend(loc='lower left')
    plt.show()

    for (lbl, fn) in funs:
        E = []
        for i in N:
            xs = linspace(0, x0_2, i+1)
            fehler = abs(fn(xs, 0, x0_2, f) - I0_2)
            E.append(fehler)
        plt.semilogy(N, E, label=lbl)

    I = romberg(0, x0_2, 4)
    plt.semilogy([4, 8, 16, 32], [abs(I[i,i] - I0_2) for i in range(4)], label='Romberg')

    plt.legend(loc='lower left')
    plt.show()
    
a7_2_5()

'''Intervalen bei Romberg

Die Eingaben bei Romberg werden benutzt um die Intervallängen zu bestimmen. dies funktionert anders
als die *Regel. Bei I[i,0] wird eine Intervallänge von x0 / 4*2**i benutzt. Also ist die Anzahl der
Intervale n = 4*2**i. Das heißt, um 40 Intervalle zu verwenden, brauchen wir log(40/4) = 3.32. Wir
setzen 4 in die Romberg-Folge, und I[i,i] ist was wir für Intervallänge n = 4*2**i verwenden.

h = x0 / 4*2**i, und h = (b - a)/n mit b = x0 und a = 0

daraus folgt 2**i = b*n / 4b -> 2**i = n/4

'''

'''
Bemerkung zu A7.2.5
    Bei x0_1= pi/2.0:
        * Trapez: bessere Genauigkeit als Simpson, aber schlechter als der Rest
        * Simpson: hat hoehere Fehler als die anderen Methoden
        * h-Extrapolation:
        * Gauss-Quadratur: hat hoechste Genauigkeit von allen Methoden, ab n = 10 ist der Fehler vernachlaessigbar klein.

    Bei x0_2 = pi
        * Trapez: hat eine hohe Genauigkeit wie bei gauss-quadratur Verfahren, viel besser als bei x0_1
        * Simpson: hat hoehere Fehler als die anderen Methoden
        * h-Extrapolation: 
        * Gauss-Quadratur: hat hoechste Genauigkeit von allen Methoden, ab n = 15 ist der Fehler vernachlaessigbar klein. 
    
    Aufgrund der Fehler, die sehr klein sind und deswegen als 0 resultieren (bei Trapez und Gauss). Da wird kein ergebnis gezeigt, 
    deswegen wird eine Zahl dementsprechend gerundet, damit es gezeigt wird.
'''
