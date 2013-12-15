#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import exp
import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt

def F(tn,y,alpha,R=1,b=0.01):
    return alpha*(y**2)*(y-R)-b*y

def euler(F,a,b,h,y0,alpha):
    yn = y0
    tn = a
    n = int(round((b-a) / h))

    x = [a]
    y = [y0]
    
    while n > 0:
        yn = yn + h*F(tn,yn,alpha)
        tn += h
        
        n -= 1
        
        x.append(tn)
        y.append(yn)
    return x,y

#Aufgabe 8.2
def plot_8_2_a():
    alphas = np.arange(0.001,1,0.05)
    a = 0
    b = 2000
    y0 = 1
    h = 0.1
    for alpha in alphas:
        ts = []
        ys = []
        ts,ys = euler(F, a, b, h, y0, alpha)
        plt.plot(ts,ys,label='alpha= '+str(alpha))
    #Aufgabe 8.2c
    ts = []
    ys = []
    a_min = -0.04
    ts,ys = euler(F, a, b, h, y0,a_min)
    plt.plot(ts,ys,'.',label='a_min = '+str(a_min))
    a_min_minus = a_min -0.01
    ts,ys = euler(F, a, b, h, y0, a_min_minus)
    plt.plot(ts,ys,label='a_min - 0.01')
    a_min_plus = a_min + 0.01
    ts,ys = euler(F, a, b, h, y0, a_min_plus)
    plt.plot(ts,ys,label='alpha + 0.01')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    plot_8_2_a()
