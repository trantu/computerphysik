#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import cos, pi, sqrt
from itertools import cycle, islice
from numpy import linspace

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

a, b, n = 0, pi, 21
xs = linspace(a, b, n)
print '### x0 = pi'
for (lbl, fn) in funs:
    print '%s: %f' % (lbl, fn(xs, a, b, f))

# Aufgabe 7.2.4    
from scipy.special.orthogonal import p_roots
def gauss_quadratur(a,b,f,n):
    
    [x,w] = p_roots(n)
    summe = 0.0
    print len(x)
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
# Aufgabe 7.2.4 
   
from math import exp,cos    

print gauss_quadratur(-1, 1, exp, 20)
