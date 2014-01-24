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
    vorfaktor = (b-a)/N
    summe = 0.0
    r = RandomGenerator()
    for i in range(1,N):
        r =  float(r.next())/(2**31 - 1)
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
            
def mittelpunkt(f,a,b,N):
    vorfaktor = np.float((b-a)/N)
    summe = 0.0
    for i in range(1,N):
        summe = summe + f(a+(b-a)*(i-0.5)/np.float(N))
    return vorfaktor*summe
def u12_2():
    N = 2**16
    ns = [2**(x) for x in range(-1,17)]
if __name__ == '__main__':
    r = RandomGenerator()
    for _ in range(10):
        print float(r.next())/(2**31 - 1)
#!
