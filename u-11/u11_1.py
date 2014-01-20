# -*- coding: utf-8 -*-
'''
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
    plt.xlabel('0.1 tau = 1 Datenpunkt')
    plt.ylabel('Phi(tau)')
    plt.show()

u11_1_1()
'''
Zeitmessung:
phi(tau = 0.0): time =  3.07 s
phi(tau = 10.0): time =  2.97 s
'''
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
    for i in range(2*l):
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
            print 'print: ', j
            printed = 1
    elapsed = (time.clock() - start)
    print 'u11.1.2 Zeitmessung fuer alle tau bis 10000: ', elapsed
    #print phis
    plt.plot(range(10000), phis[0:10000])
    plt.title('u11_1_2')
    plt.xlabel('tau in Ps (0.1 tau = 10 Datenpunkte)')
    plt.ylabel('Phi\'s')
    plt.show()
    #summe = 0.0
    #for i in range(l):
    #    summe = summe + Px[i]+Py[i]+Pz[i]
    #print summe.real
    #print np.fft.ifft(np.fft.fft(data,30,0))
    
u11_1_2()
'''
u11_1_2 
Zeitmessung:
fuer alle tau (1:1000) -> Zeit = 0.75 s
'''

def phi(tau = 10000):
    Px = np.zeros(2*l)
    Py = np.zeros(2*l)
    Pz = np.zeros(2*l)
    Px[0:l] = data[:,0]
    Py[0:l] = data[:,1]
    Pz[0:l] = data[:,2]
    Px = np.fft.fft(Px)
    Py = np.fft.fft(Py)
    Pz = np.fft.fft(Pz)
    
    for i in range(2*l):
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
    def X(f_):
        f = (10**(-4))*f_ # 1 Ghz = 10**9 1/s -> 0.1 Ghz = 10**(-4) 1/ps
        vorfactor = 1.8769
        z1 = (2*np.pi*f)*1j
        a = 0
        b = 1046 # Obergrenze abgelesen von u11_1_2()
        summe = 0.0
        for t in range(a,b+1):
            summe = summe + np.exp(-2*np.pi*f*(t/10.0)*1j)*phis[t]
        return vorfactor*(phis[0] - z1*(summe))
    ys_real = [X(freq).real for freq in np.arange(0.1,100.1,0.1)]
    print 'real 0: ', ys_real
    ys_imag = [X(freq).imag for freq in np.arange(0.1,100.1,0.1)]
    plt.figure()
    plt.title('u11.1.3 Real-Teil-Plot')
    plt.semilogy(ys_real,np.arange(0.1,100.1,0.1))
    plt.ylabel('Frequenzen')
    plt.xlabel('X_real')
    plt.figure()
    plt.title('u11.1.3 Imaginaer-Teil-Plot')
    plt.semilogy(ys_imag,np.arange(0.1,100.1,0.1))
    plt.ylabel('Frequenzen')
    plt.xlabel('X_imaginaer')
    plt.show()
u11_1_3()
