#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

from math import exp
import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt

def F(y,alpha,R=1,beta=0.01):
    return alpha*(y**2)*(y-R)-beta*y

def euler(F,a,b,h,y0,alpha):
    yn = y0
    tn = a
    n = int(round((b-a) / h))

    x = [a]
    y = [y0]
    
    #print("[EULER] tn = %s yn = %s"%(tn, yn))
    while n > 0:
        yn = yn + h*F(yn,alpha)
        tn += h
        
        n -= 1
        
        x.append(tn)
        y.append(yn)
        
        #print("[EULER] tn = %f yn = %s"%(tn, yn))
        
    return x,y

#Aufgabe 8.2
def plot_8_2_a():
    alphas = np.arange(1,10,1)
    a = 0
    b = 25
    y0 = 1
    h = 0.001
    for alpha in alphas:
        ts,ys = euler(F, a, b, h, y0, alpha)
        plt.plot(ts,ys,label='alpha= '+str(alpha))
    #Aufgabe 8.2c
    a_min = 0.002
    ts,ys = euler(F, a, b, h, y0,a_min)
    plt.plot(ts,ys,label='a_min')
    a_min_minus = a_min -0.01
    ts,ys = euler(F, a, b, h, y0, a_min_minus)
    plt.plot(ts,ys,label='a_min - 0.01')
    a_min_plus = a_min + 0.01
    ts,ys = euler(F, a, b, h, y0, a_min_plus)
    plt.plot(ts,ys,label='alpha + 0.01')
    plt.legend()
    plt.show()

def rungeKutta(F,a,b,h,y0):
    yn = y0
    tn = a
    n = int(round((b-a) / h))
    #print("[RUNGE KUTTA] tn = %f yn = %s"%(tn, yn))
    
    x = [a]
    y = [y0]
    
    
    while n > 0:
        k1 = F(tn,yn)
        k2 = F(tn + 0.5*h , yn + h*0.5*k1)
        k3 = F(tn + 0.5*h , yn + h*0.5*k2)
        k4 = F(tn + h , yn + h*k3)
        
        yn = yn + h/6.0*(k1 + 2*k2 + 2*k3 + k4)
        tn += h
        
        n -= 1
        
        # Ergebnisse speichern
        x.append(tn)
        y.append(yn)
        #print("[RUNGE KUTTA] tn = %f yn = %s"%(tn, yn))
        
    return [x,y]
if __name__ == '__main__':
    plot_8_2_a()
