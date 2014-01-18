# -*- coding: utf-8 -*-
'''
Created on Nov 29, 2013

@author: Carlos Matin Nieto, Tu Tran
'''
import time
import numpy as np
import matplotlib.pyplot as plt
# Daten einlesen
data = np.recfromtxt('P.txt')
(l,_) = np.shape(data)
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
    plt.title('u11_1_1')
    plt.show()

u11_1_1()

#u11_1_2
def u11_1_2():
    start = time.clock()
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
    tau = 10000
    phis = []
    printed = 0
    for j in range(tau):
        summe = 0.0
        #print range(l-j)
        summe = summe + Px[j]+Py[j]+Pz[j]
        normierung = summe.real*(l-j)**(-1)
        phis.append(normierung)
        #print normierung     
        # erste Nullstelle fuer u11_1_3 finden (Obergrenze des Intergrals)
        if normierung >= 0.0 and normierung <= 0.001 and not printed:
            print 'print', j
            printed = 1
    elapsed = (time.clock() - start)
    print 'Time u11.1.2: ', elapsed
    #print phis
    plt.plot(range(1000), phis[0:1000])
    plt.title('u11_1_2')
    plt.show()
    # TODO normieren, Zeit messen
    #summe = 0.0
    #for i in range(l):
    #    summe = summe + Px[i]+Py[i]+Pz[i]
    #print summe.real
    #print np.fft.ifft(np.fft.fft(data,30,0))
    
u11_1_2()

def phi(tau = 100000):
    Px = np.zeros(2*l)
    Py = np.zeros(2*l)
    Pz = np.zeros(2*l)
    Px[0:l] = data[:,0]
    Py[0:l] = data[:,1]
    Pz[0:l] = data[:,2]
    Px = np.fft.fft(Px)
    Py = np.fft.fft(Py)
    Pz = np.fft.fft(Pz)
    
    for i in range(l):
        Px[i] = abs(Px[i]*Px[i])
        Py[i] = abs(Py[i]*Py[i])
        Pz[i] = abs(Pz[i]*Pz[i])
        
    Px = np.fft.ifft(Px)
    Py = np.fft.ifft(Py)
    Pz = np.fft.ifft(Pz)
    phis = []
    for j in range(tau):
        summe = 0.0
        #print range(l-j)
        summe = summe + Px[j]+Py[j]+Pz[j]
        normierung = summe.real*(l-j)**(-1)
        phis.append(normierung)
    return phis
def u11_1_3():
    phis = phi()
    def X(f):
        vorfactor = 1.8769
        phi_0 = phis[0]
        z1 = -(2*np.pi*f)*1j
        a = 0
        b = 1046 # Obergrenze abgelesen von u11_1_2()
        summe = 0.0
        for t in range(a,b):
            summe = summe + np.exp(-2*np.pi*f*t*1j)*phis[t]
        
        return vorfactor*(phi_0 - z1*(summe))
    ys_real = [X(freq).real for freq in np.arange(0.1,100,0.1)]
    ys_imag = [X(freq).imag for freq in np.arange(0.1,100,0.1)]
    plt.figure()
    plt.title('u11.1.3 Real-Teil-Plot')
    plt.semilogy(np.arange(0.1,100,0.1),ys_real)
    plt.figure()
    plt.title('u11.1.3 Imaginaer-Teil-Plot')
    plt.semilogy(np.arange(0.1,100,0.1),ys_imag)
    plt.show()
u11_1_3()
