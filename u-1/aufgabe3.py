#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto
# Tu Tran

import time
import math
import numpy
import matplotlib.pyplot as plot

# cos through Taylor
def tcos(x, N):
    fn = lambda x, n: pow(-1, n) * (pow(x, 2*n) / math.factorial(2*n))
    return math.fsum([fn(x, n) for n in range(0, N)])

# Return a tuple with the result of the function and how long it took
def bench(fn):
    start = time.clock()
    res = fn()
    stop = time.clock()

    return (res, stop - start)

## First some sanity checking to see that our numbers look rougnly the same
res = [(numpy.cos(n), tcos(n, 15)) for n in [0.0, 1.0, math.pi, 2*math.pi]]

for val in res:
       print("%f\t%f" % val)

## Now let's compare the results of tcos with numpy.cos

# These lists are constant for all of the plots
xs = numpy.arange(0, 4*math.pi, 0.1)
theirs = [numpy.cos(n) for n in xs]

# Return tcos(x, N) values for our list
def splot(n):
    return [tcos(x, n) for x in xs]

# Calculate the pairweise difference between two lists
def diff(a, b):
    return [abs(a - b) for (a, b) in zip(ours, theirs)]

# Store the instructions for the plots. Each tuple stores N and a
# tuple with the subplot we want to draw it on
plots = [(1, (2,3,1)), (3, (2,3,2)), (5, (2,3,3)), (7, (2,3,4)),
         (9, (2,3,5)), (11, (2,3,6))]
for (n, location) in plots:
    plot.subplot(*location)
    plot.title("n = %d" % n)

    ours = splot(n)
    plot.plot(xs, ours)
    plot.plot(xs, theirs)
    plot.plot(xs, diff(ours, theirs))

    plot.xlabel(u'x')
    plot.ylabel(u'cos(x)')

    plot.legend((u'$t\cos_N(x)$', u'$\cos(x)$', u'Fehler'), loc='lower left')
    plot.axis([0, 4*math.pi, -1, 1])

plot.show()

## Let's benchmark it

rng = range(0, 1000, 10)
ours = [bench(lambda: tcos(0, n)) for n in rng]
theirs = [bench(lambda: numpy.cos(n)) for n in rng]
plot.subplot(111)
plot.plot(rng, [t for (_, t) in ours])
plot.plot(rng, [t for (_, t) in theirs])
plot.show()
