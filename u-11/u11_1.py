# -*- coding: utf-8 -*-
'''
Created on Nov 29, 2013

@author: Carlos Matin Nieto, Tu Tran
'''
import time
import numpy as np
import matplotlib.pyplot as plt

# Aufgabe 11.1.1
def autokorrelation(T, tao, Pdata):
    autokorr = 0.0
    tao_ = int(tao * 10)
    print tao_
    schritte = T - tao_ # 0.1 ps ist eine Zeile
    t = 0
    normierung = (T - tao_)**(-1) 
    while t + tao_ < schritte:
        #print int(tao_ + t)
        autokorr = autokorr + Pdata[t][0]*Pdata[t+tao_][0]+Pdata[t][1]*Pdata[t+tao_][1]+Pdata[t][2]*Pdata[t+tao_][2]
        t = t + 1
        #print t
        #print autokorr
    
    return normierung*autokorr

def u11_1_1():
    data = np.recfromtxt('P.txt')
    (l,_) = np.shape(data)
    start = time.clock()
    print autokorrelation(l, 0.0, data)
    elapsed = (time.clock() - start)
    print 'tau = 0.0 time = ',elapsed
    start = time.clock()
    print autokorrelation(l, 10.0, data)
    elapsed = (time.clock() - start)
    print 'tau = 10.0 time = ',elapsed
    phis = []
    for i in range(0,100,5):
        phis.append(autokorrelation(l, i, data))
    plt.plot(range(0,100,5),phis)
    plt.show()
# 11.1.1
u11_1_1()
def u11_1_2():
    data = np.recfromtxt('P.txt')
    (l,_) = np.shape(data)
    Px = np.zeros(2*l)
    Py = np.zeros(2*l)
    Pz = np.zeros(2*l)
    Px[0:l] = data[:,0]
    Py[0:l] = data[:,1]
    Pz[0:l] = data[:,2]
    Px = np.fft.fft(Px)
    Py = np.fft.fft(Py)
    Pz = np.fft.fft(Pz)
    #print z[99997]
    for i in range(l):
        Px[i] = abs(Px[i]*Px[i])
        Py[i] = abs(Py[i]*Py[i])
        Pz[i] = abs(Pz[i]*Pz[i])
        
    Px = np.fft.ifft(Px)
    Py = np.fft.ifft(Py)
    Pz = np.fft.ifft(Pz)
    print np.shape(Px)
    print np.shape(Py)
    print np.shape(Pz)
    tao = 10000
    phis = []
    for j in range(tao):
        summe = 0.0
        #print range(l-j)
        summe = summe + Px[j]+Py[j]+Pz[j]
        phis.append(summe.real)
    print 'bla'
    plt.plot(range(1000), phis[0:1000])
    plt.show()
    #summe = 0.0
    #for i in range(l):
    #    summe = summe + Px[i]+Py[i]+Pz[i]
    #print summe.real
    #print np.fft.ifft(np.fft.fft(data,30,0))
    
u11_1_2()
