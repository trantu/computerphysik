#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

import math
import numpy
import matplotlib.pyplot as plot

def cos_1_ableitung(x):
    return -1*math.sin(x);
def cos_2_ableitung(x):
    return -1*math.cos(x);

def diff_quot(fn, xs, ys):
    values = [(ys[i+1]-ys[i-1])/(xs[i+1]-xs[i-1]) for i in range(1, len(xs) -1)]
    control = map(fn, xs[1:-1])
    results = zip(values, control)
    errors = [abs(a - b) for (a, b) in results]
    return (results, errors)

def delta_X(k):
    return 2*math.pi/2**k
    
def k(f, f_ableitung):
    rg = range(10)
    rg.reverse()
    delta_x = map(delta_X,rg)
    #print delta_x
    counter = len(rg)-1 
    
    for i in delta_x:
        #print i
        xs = numpy.arange(0.0,20+i,i)
        ys = map(f, xs)
        (res,abs_fails) = diff_quot(f_ableitung, xs, ys) 
        #print abs_fails
        for j in abs_fails:
            print "fehler = %g" %j
            if j > 0.01:
                print "Fehler = %g , Delta_x = %g" %(j,i)
                return counter
        counter -= 1           
    return counter
print "k = %i" % k(math.cos, cos_1_ableitung)

'''
Ausgabe:
Ab k = 4 mit Delta_x = 0.392699 wird der Fehler bei 0.0180345 groesser sein als 0.01. 
Wenn k >= 5 wird der Fehler kleiner sein als 0.01, deswegen waehlen wir k = 5 fuer delta_x = 0.196349540849
'''
# k = 5
delta_x_k = 2*math.pi/2**5
#Intervall fuer normale Funktion
interval = numpy.arange(-math.pi, math.pi+delta_x_k, delta_x_k)
#Intervall fuer diff_quotient der 1ten Ableitung
interval_diff = interval[1:-1]
#Intervall fuer diff_quotient der 2ten Ableitung
interval_diff_2te_Ableitung = interval_diff[1:-1]
# 1. Ableitung von cos soll durch diff_quotient angenaehert werden
ys = map(math.cos, interval)
res1, abs_fehler = diff_quot(cos_1_ableitung, interval, ys)

res_a1, res_b1 = zip(*res1)

plot.plot(interval_diff, res_b1)
plot.plot(interval_diff, res_a1)

# 2. Ableitung von cos soll durch diff_quotient angenaehert werden
ys = res_a1
res2, abs_fehler2 = diff_quot(cos_2_ableitung, interval_diff, ys)

res_a2, res_b2 = zip(*res2)

plot.plot(interval_diff_2te_Ableitung, res_a2)
plot.plot(interval_diff_2te_Ableitung, res_b2)

plot.legend(("$\cos'(x)$", "$-\sin(x)$", "$-\sin'(x)$", "$-\cos(x)$"))
plot.xlabel("$x$")
plot.ylabel("$y$")

plot.show()

