#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import cos, pi, sqrt, log10
from itertools import cycle, islice
from numpy import linspace
import matplotlib.pyplot as plt
from scipy.special.orthogonal import p_roots
import numpy as np
from numpy.core.numeric import dtype

def take(seq, n):
    return list(islice(seq, 0, n))

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
    s = sum([f(xs[i]) for i in range(1, n)])

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
    ('Gauss', gauss_quadratur),
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
#print 'Gauss-Quadratur', gauss_quadratur(a, b, f, n)
#Aufgabe 7.2.5
def a7_2_5():
    N = range(2,41,2)
    I0_1 = round(pi/(3*sqrt(3)),14)
    I0_2 = round(pi/sqrt(3),13)
    x0_1= pi/2.0
    x0_2 = pi
    for (lbl, fn) in funs:
        if lbl != 'Genau' and lbl != 'Rechteck':
            E = []
            for i in N:
                xs = linspace(0, x0_1, i)
                fehler=abs(fn(xs, 0, x0_1, f)-I0_1)
                E.append(fehler)
            plt.semilogy(N,E,label=lbl)
            plt.legend()
    plt.show()
    print abs(log10(0.5))
    for (lbl, fn) in funs:
        if lbl != 'Genau' and lbl != 'Rechteck':
            E = []
            for i in N:
                xs = linspace(0, x0_2, i)
                fehler = abs(fn(xs, 0, x0_2, f)-I0_2)
                '''
                if fehler == 0.0:
                    print i
                    print lbl
                    print I0_2
                    print fn(xs, 0, x0_2, f)
                '''
                E.append(fehler)
            #print E
            plt.semilogy(N,E,label=lbl)
            plt.legend()
    plt.show()
    
a7_2_5()

'''
Bemerkung zu A7.2.5
    Bei x0_1= pi/2.0:
        * Trapez: bessere Genauigkeit als Simpson, aber schlechter als der Rest
        * Simpson: hat hoehere Fehler als die anderen Methoden
        * h-Extrapolation:
        * Gauss-Quadratur: hat hoechste Genauigkeit von allen Methoden, ab n = 10 ist der Fehler vernachlaessigbar klein.
    Aufgrund der Fehler, die sehr klein sind und deswegen als 0 resultieren
    Bei x0_2 = pi
        * Trapez: hat eine hohe Genauigkeit wie bei gauss-quadratur Verfahren, viel besser als bei x0_1
        * Simpson: hat hoehere Fehler als die anderen Methoden
        * h-Extrapolation: 
        * Gauss-Quadratur: hat hoechste Genauigkeit von allen Methoden, ab n = 15 ist der Fehler vernachlaessigbar klein. 
'''
