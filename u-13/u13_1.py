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
def fg(y):
    return 0.5 * y/(np.sqrt(np.log(1+(y**2)* 0.25)))
def u13_1_1():
    def run(N):
        delta = 10**(-3)
        ergebnisse = []
        b = 1
        R = 0.0
        zero = 10**-8 # weil x = 0 undefiniert ist
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
    Bei 1000 Schüssen:
    R = 0.0515669030678
    0.0715034424513
    0.081471712143
    0.0715034424513
    0.0715034424513
    0.0415986333761
    0.0515669030678
    0.0316303636844
    0.0515669030678
    0.0515669030678
    [2.1634620301917322, 2.122600313525389, 2.2237659672917136, 2.122600313525389, 2.1833985695751368, 1.7887042242015414, 1.76827336586837, 2.2347228748829493, 2.406655054390725, 2.1026637741419845]
    Bei 100000 Schüssen:
    R = 0.0589434226397
    0.0620335862441
    0.0624323170318
    0.060937076578
    0.0653231152424
    0.0652234325455
    0.0607377111842
    0.0627313651226
    0.0607377111842
    0.062232951638
    [2.1316236746115047, 2.1626810359988169, 2.1968227988940945, 2.1442570233585512, 2.1398273148956868, 2.1518872834087195, 2.1586492394166563, 2.1667227189599716, 2.1534813876524277, 2.1397770640937464]
    '''
def u13_1_2():
    def run(N):
        ergebnisse = []
        b = 2
        zero = 10**-6 # weil y = 0 undefiniert ist
        for _ in range(10):
            erg = hit_miss_monte(fg, fg(zero), fg(b),N, zero, b)
            ergebnisse.append(erg)
        print ergebnisse
    run(1000)
    run(100000)
    '''
    Kommentar:
    Bei 1000 Schüssen:
    [2.1330824975505389, 2.143945502349192, 2.1387151667053961, 2.1515898390593549, 2.1479688374598038, 2.1531991731035998, 2.1511875055482936, 2.1507851720372329, 2.1435431688381303, 2.141129167771763]
    Bei 100000 Schüssen:
    [2.1447380993659828, 2.1448266127384157, 2.1435270754976883, 2.145438159675229, 2.1448064960628632, 2.1454462063454502, 2.1436115655350112, 2.1440742490727316, 2.1453053896165786, 2.144814542733084]
    
    Die Werte hier im Vergleich zu u13.1.1 sind sehr stabil, besonders bei 10**5 Schüssen. Hier sieht man, dass es sehr schnell konvergiert.
    In u13.1.1 selbst bei 10**5 Schüssen haben die Werte immer noch hohe Schwankungen.
    '''   
if __name__ == '__main__':
    u13_1_1()
    u13_1_2()
