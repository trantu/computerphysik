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

def rungeKutta(F1,F2,F3,F4,a,b,h,y0):
    yn = y0
    tn = a
    n = int(round((b-a) / h))
    #print("[RUNGE KUTTA] tn = %f yn = %s"%(tn, yn))
    
    x = [a]
    y = [y0]
    
    
    while n > 0:
        k1_1 = F1(tn,yn[0],yn[1],yn[2],yn[3])
        k1_2 = F2(tn,yn[0],yn[1],yn[2],yn[3])
        k1_3 = F3(tn,yn[0],yn[1],yn[2],yn[3])
        k1_4 = F4(tn,yn[0],yn[1],yn[2],yn[3])
        
        k2_1 = F(tn + 0.5*h,yn[0] + h*0.5*k1_1,yn[1] + h*0.5*k1_1,yn[2] + h*0.5*k1_1,yn[3] + h*0.5*k1_1)
        k2_2 = F(tn + 0.5*h,yn[0] + h*0.5*k1_2,yn[1] + h*0.5*k1_2,yn[2] + h*0.5*k1_2,yn[3] + h*0.5*k1_2)
        k2_3 = F(tn + 0.5*h,yn[0] + h*0.5*k1_3,yn[1] + h*0.5*k1_3,yn[2] + h*0.5*k1_3,yn[3] + h*0.5*k1_3)
        k2_4 = F(tn + 0.5*h,yn[0] + h*0.5*k1_4,yn[1] + h*0.5*k1_4,yn[2] + h*0.5*k1_4,yn[3] + h*0.5*k1_4)
        
        k3_1 = F(tn + 0.5*h,yn[0] + h*0.5*k2_1,yn[1] + h*0.5*k2_1,yn[2] + h*0.5*k2_1,yn[3] + h*0.5*k2_1)
        k3_2 = F(tn + 0.5*h,yn[0] + h*0.5*k2_2,yn[1] + h*0.5*k2_2,yn[2] + h*0.5*k2_2,yn[3] + h*0.5*k2_2)
        k3_3 = F(tn + 0.5*h,yn[0] + h*0.5*k2_3,yn[1] + h*0.5*k2_3,yn[2] + h*0.5*k2_3,yn[3] + h*0.5*k2_3)
        k3_4 = F(tn + 0.5*h,yn[0] + h*0.5*k2_4,yn[1] + h*0.5*k2_4,yn[2] + h*0.5*k2_4,yn[3] + h*0.5*k2_4)
        
        k4_1 = F(tn + h,yn[0] + h*k3_1,yn[1] + h*k3_1,yn[2] + h*k3_1,yn[3] + h*k3_1)
        k4_2 = F(tn + h,yn[0] + h*k3_2,yn[1] + h*k3_2,yn[2] + h*k3_2,yn[3] + h*k3_2)
        k4_3 = F(tn + h,yn[0] + h*k3_3,yn[1] + h*k3_3,yn[2] + h*k3_3,yn[3] + h*k3_3)
        k4_4 = F(tn + h,yn[0] + h*k3_4,yn[1] + h*k3_4,yn[2] + h*k3_4,yn[3] + h*k3_4)

        
        yn[0] = yn[0] + h/6.0*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1)
        yn[1] = yn[1] + h/6.0*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)
        yn[2] = yn[2] + h/6.0*(k1_3 + 2*k2_3 + 2*k3_3 + k4_3)
        yn[3] = yn[3] + h/6.0*(k1_4 + 2*k2_4 + 2*k3_4 + k4_4)
        tn += h
        
        n -= 1
        
        # Ergebnisse speichern
        x.append(tn)
        y.append(yn)
        #print("[RUNGE KUTTA] tn = %f yn = %s"%(tn, yn))
        
    return [x,y]
if __name__ == '__main__':
    plot_8_2_a()
