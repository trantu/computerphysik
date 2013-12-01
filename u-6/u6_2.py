'''
Created on Nov 29, 2013

@author: Carlos Matin Nieto, Tu Tran
'''
# -*- coding: utf-8 -*-
import numpy as np
from numpy import * 
from copy import deepcopy
import matplotlib.pyplot as plt
import operator as op

def f(a0,a1,a2,a3,t):
    return (exp(-a0 * t )) * (a2*sin(a1*t) + a3*cos(a1*t))
def g(a0,a1,a2,t):
    return (exp(-a0 * t)) * a2 * sin(a1 * t)

#def gauss_newton(jacob):
def jacob_f(xs,aS):
    
    def f_a0(x):
        return (-x*exp(-aS[0]*x))* (aS[2]*sin(aS[1]*x) + aS[3]*cos(aS[1]*x))         
    def f_a1(x):
        return (x*exp(-aS[0]*x))* (aS[2]*cos(aS[1]*x) - aS[3] * sin(aS[1]*x))
    def f_a2(x):
        return (exp(-aS[0]*x))* sin(aS[1]*x)
    def f_a3(x):
        return (exp(-aS[0]*x))* cos(aS[1]*x)
    jmatrix_f = np.matrix(np.zeros((len(xs), 4)))
    for i in range(len(xs)):
        jmatrix_f[i] = [f_a0(xs[i]), f_a1(xs[i]), f_a2(xs[i]), f_a3(xs[i])]
        
    return jmatrix_f

def jacob_g(xs,aS):
    
    def g_a0(x):
        return (-x*exp(-aS[0]*x))* aS[2]*sin(aS[1]*x)         
    def g_a1(x):
        return (x*exp(-aS[0]*x))* aS[2]*cos(aS[1]*x)
    def g_a2(x):
        return (exp(-aS[0]*x))*sin(aS[1]*x)
    
    jmatrix_g = np.matrix(np.zeros((len(xs), 3)))
    for i in range(len(xs)):
        jmatrix_g[i] = [g_a0(xs[i]), g_a1(xs[i]), g_a2(xs[i])]
        
    return jmatrix_g

def g_f(xs,ys,a_s):
    #print np.exp(-a_s[0] * xs )
    return (ys - (np.exp(-a_s[0] * xs )) * (a_s[2]*np.sin(a_s[1]*xs) + a_s[3]*np.cos(a_s[1]*xs)))

def g_g(xs,ys,a_s):
    return (ys - (np.exp(-a_s[0] * xs)) * a_s[2] * np.sin(a_s[1] * xs))

def gauss_newton(aS,g,jacobi):
    _,spalten = np.shape(jacobi)
    #print aS
    #if((np.shape(g) != np.shape(jacobi)) or len(aS.getT()) != spalten):
     #   raise ValueError('dimension is not equal')
    #n = 1
    #a = aS
    #while(n < 101):
        
    delta = np.linalg.solve(jacobi.getT() * jacobi , jacobi.getT() * np.matrix(g).getT())
    delta = delta.getT().tolist()[0]
    
    
    #print("Schritt %i delta = %s "%(n, delta))
    
    aS = aS + delta
    #print delta
        #n += 1
    return aS,delta


def gauss_newton_modi(aS,g,jacobi,f,xs,ys):
    _,spalten = np.shape(jacobi)
        
    delta = np.linalg.solve(jacobi.getT() * jacobi , jacobi.getT() * np.matrix(g).getT())
    delta = delta.getT().tolist()[0]
    np.linalg.norm(delta,2)
    return find_min_k(aS, g, f, np.array(delta), xs, ys)

def find_min_k(a_n,g_n,f,delta,xs,ys):
    k = 0
    #print delta
    while(1):
        a_n1 = a_n + (2**(-k))* delta
        g_n1 = f(xs,ys,a_n1)
        E = np.linalg.norm(g_n, 2)**2
        E1 = np.linalg.norm(g_n1,2)**2
        if E > E1 or 2**(-k) < 10**(-6):
            break
        k += 1 
    return k,a_n1,delta
    
def main():
    a1 = [0.8, 6.4, 4.2, -0.3]
    a2 = [0.3, 5.4, 7.2, -1.3]
    a3 = [1.0, 7.0, -6.0, 3.0]
    b1 = [0.8, 6.4, 4.2]
    b2 = [0.3, 5.4, 7.2]
    b3 = [1.0, 7.0, -6.0]
    
    
    data = np.recfromtxt('data2.txt')
    xs = data[:,0]
    ys = data[:,1]
    xs_erstes_viertel = xs[:len(xs)/4]
    
    def run_gauss_newton_g(aS, xs, ys):
        a = np.array(aS)
        #print 'as: ', a
        c = 0
        while c < 200:
            gg = g_g((np.array(xs)),(np.array(ys)),a)
            #print np.linalg.norm(gg, 2)
            jacobi = jacob_g(xs, a)
            a,delta = gauss_newton(a, gg, jacobi)
            #print 'delta: ', delta
            #print 'a: ', a 
            #sa = map(op.add, delta, a)
            if(np.linalg.norm(delta,2) < (10**(-6))):
                break
            c +=1
        print 'Iterationsschritte: ', c
        return a   
    def run_gauss_newton_f(aS, xs, ys):
        a = np.array(aS)
        #print 'as: ', a
        c = 0
        while c < 200:
            gf = g_f((np.array(xs)),(np.array(ys)),a)
            #print np.linalg.norm(gf, 2)
            jacobi = jacob_f(xs, a)
            a,delta = gauss_newton(a, gf, jacobi)
            #print 'f::: ',a
            #print 'delta: ', delta
            #print 'a: ', a 
            #sa = map(op.add, delta, a)
            if(np.linalg.norm(delta,2) < (10**(-6))):
                break
            c +=1
        print 'Iterationsschritte: ', c
        return a 
    def run_gauss_newton_modi_f(aS, xs, ys):
        a = np.array(aS)
        #print 'as: ', a
        c = 0

        while c < 100:
            gf = g_f((np.array(xs)),(np.array(ys)),a)
            #print np.linalg.norm(gf, 2)
            jacobi = jacob_f(xs, a)
            k,a,delta = gauss_newton_modi(a, gf, jacobi, g_f, np.array(xs), np.array(ys))
            if(np.linalg.norm(delta,2) < (10**(-6))):
                break                                            
            #print 'delta: ', delta
            #sa = map(op.add, delta, a)
            #if(np.linalg.norm(delta,2) < (10**(-6))):
             #   break
            c +=1
        print 'Iterationsschritte: ', c
        return a
    
    def run_gauss_newton_modi_g(aS, xs, ys):
        a = np.array(aS)
        #print 'as: ', a
        c = 0

        while c < 100:
            gg = g_g((np.array(xs)),(np.array(ys)),a)
            #print np.linalg.norm(gf, 2)
            jacobi = jacob_g(xs, a)
            k,a,delta = gauss_newton_modi(a, gg, jacobi, g_g, np.array(xs), np.array(ys))
            if (np.linalg.norm((2**(-k))* delta,2) < 10**(-6)): break                                            
            #print 'delta: ', delta
            #sa = map(op.add, delta, a)
            #if(np.linalg.norm(delta,2) < (10**(-6))):
             #   break
            c +=1
        print 'Iterationsschritte: ', c
        return a
 
    # Running

    print u'Ungedaempten Verfahren'
    as_ = run_gauss_newton_f(a1, xs, ys)
    print 'f an a1: ', as_
    print 'f an a2: ', run_gauss_newton_modi_f(a2, xs, ys)
    print 'f an a3: ',run_gauss_newton_f(a3, xs, ys)
    print 'g an b1: ',run_gauss_newton_g(b1, xs, ys)
    print 'g an b2: ',run_gauss_newton_g(b2, xs, ys)
    print 'g an b3: ',run_gauss_newton_g(b3, xs, ys)

    print u'Gedaempten Verfahren'
    print 'f an a1: ',run_gauss_newton_modi_f(a1, xs, ys)
    print 'f an a2: ',run_gauss_newton_modi_f(a2, xs, ys)
    print 'f an a3: ',run_gauss_newton_modi_f(a3, xs, ys)
    print 'g an b1: ',run_gauss_newton_modi_g(b1, xs, ys)
    print 'g an b2: ',run_gauss_newton_modi_g(b2, xs, ys)
    print 'g an b3: ',run_gauss_newton_modi_g(b3, xs, ys)

    plt.plot(xs, ys, label='Daten')
    pairs = [
        (xs[:len(xs)/4], ys[:len(ys)/4]),
        (xs[3*len(xs)/4:], ys[3*len(ys)/4:]),
        (xs[::5],  ys[::5]),
        (xs[::20], ys[::20]),
        (xs[::40], ys[::40])
        ]

    for (xss, yss) in pairs:
        as_ = run_gauss_newton_f(a1, xss, yss)
        ff = lambda as_, t: f(as_[0], as_[1], as_[2], as_[3], t)
        plt.plot(xs, [ff(as_, t) for t in xs])

        plt.labels()

    
if __name__ == '__main__':
    main()
