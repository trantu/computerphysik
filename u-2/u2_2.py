#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Tu Tran, Carlos Martín Nieto
'''
import math
import numpy
import matplotlib.pyplot as plot

def cos_1_ableitung(x):
    return -1*math.sin(x);
def cos_2_ableitung(x):
    return -1*math.cos(x);

def diff_quot(fn, xs, ys):
    """
    Perform differential quotients calculation, return the calculated values and the control zipped.
    """
    values = [(ys[i+1]-ys[i-1])/(xs[i+1]-xs[i-1]) for i in range(1, len(xs) - 1)]
    control = map(fn, xs)

    return zip(values, control)

def delta_X(k):
    return 2*math.pi/2**k

def k(f, f_ableitung):
    delta_x = map(delta_X,range(100))
    print delta_x
    counter = 0 
    for i in delta_x:
        #print i
        xs = numpy.arange(0.0,100+i,i)
        ys = map(f, xs)
        abs_fails = [abs(a-b) for (a,b) in diff_quot(f_ableitung,xs, ys)]
        print abs_fails
        for j in abs_fails:
            print "j = %g" %j
            if j > 0.01:
                print "Fehler = %g , Delta_x = %g" %(j,i)
                return counter
        counter += 1           
    return counter
#print res[counter]
print "k = %i" % k(math.cos,cos_1_ableitung)

'''
Fehler = 0.36338 , Delta_x = 1.5708
k = 2
'''
# k = 3 => weil damit kleiner als delta_x

delta_x_k = 2*math.pi/2**4

interval = numpy.arange(-math.pi,math.pi+delta_x_k,delta_x_k)
ys = map(math.cos, interval)
res1 = diff_quot(math.sin,interval, ys)

plot.subplot(121)
plot.plot(interval, ys)

#print res1
ys = map(lambda x: math.sin(x)*(-1), interval)
res2 = diff_quot(math.cos,interval, ys)
#print res2
resa = [a for (a,b) in res1]
#print resa
#print ys
plot.subplot(122)
plot.plot(interval, ys)
plot.show()
