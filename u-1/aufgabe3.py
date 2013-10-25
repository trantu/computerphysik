#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto
# Tu Tran

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
       print "%f\t%f" % val

## Now let's compare the results of tcos with numpy.cos

# These are always the same
xs = numpy.arange(0, 2*math.pi, 0.1)
theirs = [numpy.cos(n) for n in xs]
def splot(n):
    return [tcos(x, n) for x in xs]

plots = [(1, 231), (3, 232), (5, 233), (7, 234), (9, 235), (11, 236)]
for (n, location) in plots:
    plot.subplot(location)
    plot.title("n = %d" % n)

    ours = splot(n)
    plot.plot(xs, ours)
    plot.plot(xs, theirs)
    plot.plot(xs, [abs(a - b) for (a, b) in zip(ours, theirs)])

    plot.xlabel(u'x')
    plot.ylabel(u'cos(x)')

    plot.legend((u'$t\cos_N(x)$', u'$\cos(x)$', u'Fehler'), loc='lower left')
    plot.axis([0, 2*math.pi, -1, 1])

plot.show()
