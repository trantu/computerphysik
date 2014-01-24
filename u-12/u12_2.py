# -*- coding: utf-8 -*-
'''
Created on Nov 29, 2013

@author: Carlos Matin Nieto, Tu Tran
'''
import time,os
import numpy as np
import matplotlib.pyplot as plt

class RandomGenerator:
    # Seed the random number generator
    def __init__(self, t=time.time()):
        # Start with the magic numbers
        pid = os.getpid()
        self.a = 7**5
        self.c = 0
        self.m = 2**31 - 1

        # And seed the RNG as per formula (2)
        num = abs((64979 * t * (pid - 83)) % 104729)
        num = num / 104729
        num = num * self.m
        self.num = int(num)

    # This makes this class an iterator itself
    def __iter__(self):
        return self

    # Sensible name and lets us make this an iterator
    def next(self):
        # Calculate the new number, as per (1)
        self.num = (self.a * self.num + self.c) % self.m

        return self.num

def stuetz_monteC_integr(f,N,a,b):
    vorfaktor = (b-a)/np.float(N)
    summe = 0.0
    rs = RandomGenerator()
    for _ in range(1,N):
        r =  float(rs.next())/(2**31 - 1)
        #print summe
        summe = summe + f(a + (b-a) * r)
    return vorfaktor*summe

def hit_miss_monte(f,fmin,fmax,N,a,b):
    N_plus = 0.0
    r = RandomGenerator()
    for _ in range(N):
        ri = float(r.next())/(2**31 - 1)
        si = float(r.next())/(2**31 - 1)
        if f(a+(b-a)*ri) - fmin > (fmax-fmin)*si:
            N_plus = N_plus + 1
    return (b-a)*((N_plus/np.float(N))*(fmax - fmin) + fmin)

def f(x):
    return np.abs(np.pi*x)
def f_analytisch(n):
    return 2*n/np.pi           
def mittelpunkt(f,a,b,N):
    vorfaktor = (b-a)/np.float(N)
    summe = 0.0
    for i in range(1,N):
        summe = summe + f(a+(b-a)*(i-0.5)/np.float(N))
    return vorfaktor*summe
def u12_2():
    N = 2**16
    ns = [2**(x) for x in range(-1,17)]
    ys_stutz = [stuetz_monteC_integr(f, N, 0, n) for n in ns]
    ys_hit = [hit_miss_monte(f, 0, 1, N, 0, n) for n in ns]
    ys_mittel = [mittelpunkt(f, 0, n, N) for n in ns]
    ys_analytisch = [f_analytisch(n) for n in ns]
    
    ys_stutz_fehler = []
    ys_hit_fehler = []
    ys_mittel_fehler = []
    for i in range(len(ns)):
        ys_stutz_fehler.append((ys_stutz[i]-ys_analytisch[i])/ys_analytisch[i])
        ys_hit_fehler.append((ys_hit[i]-ys_analytisch[i])/ys_analytisch[i])
        ys_mittel_fehler.append((ys_mittel[i] - ys_analytisch[i])/ys_analytisch[i])
    print 'ys_stutz',ys_stutz
    print 'ys_hit',ys_hit
    print 'ys_mittel',ys_mittel
    print 'ys_analytisch',ys_analytisch
    print 'ys_fehler',ys_hit_fehler
    #plt.plot(ns,ys_stutz_fehler,'-*',label = 'Relativer Fehler bei Stuetz-Monte-Integration')
    plt.loglog(ns,ys_hit_fehler,'-*',label = 'Relativer Fehler Hit-miss-Monte-Integration')
    #plt.plot(ns,ys_mittel_fehler,'-*',label = 'Relativer Fehler Mittelpunkt-Integration')
    plt.show()
if __name__ == '__main__':
    u12_2()
#!
