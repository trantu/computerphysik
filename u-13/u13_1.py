# -*- coding: utf-8 -*-
'''
@author: Carlos Matin Nieto, Tu Tran
'''
import time,os
import numpy as np
import matplotlib.pyplot as plt

def hit_miss_monte(f,fmin,fmax,N,a,b):
    N_plus = 0.0
    for _ in range(N):
        ri = np.random.uniform(0,1)
        si = np.random.uniform(0,1)
        if f(a+(b-a)*ri) - fmin > (fmax-fmin)*si:
            N_plus = N_plus + 1
    return (b-a)*((N_plus/np.float(N))*(fmax - fmin) + fmin)

          
def mittelpunkt(f,a,b,N):
    vorfaktor = (b-a)/np.float(N)
    summe = 0.0
    for i in range(1,N):
        summe = summe + f(a+(b-a)*(i-0.5)/np.float(N))
    return vorfaktor*summe

def f(x):
    return 1.0/(np.sqrt(np.log(1+x)))
def u13_1_1():
    def run(N):
        delta = 10**(-3)
        ergebnisse = []
        b = 1
        R = 0.0
        zero = 10**-8 # 
        for _ in range(10):
            erg = hit_miss_monte(f, f(1), f(delta),N, delta, b)
            R = hit_miss_monte(f, f(delta), f(zero),N, zero,delta)
            ergebnisse.append(erg+R)
            print R
        print ergebnisse
    run(1000)
    run(100000)
    '''
    R hat Schwankungen bei 1000 Schüssen. Unter 100000 Schüssen sind die Werte stabiler.
    '''
    
if __name__ == '__main__':
    u13_1_1()
