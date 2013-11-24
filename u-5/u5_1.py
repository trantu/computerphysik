#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

import time
import math
import copy
import numpy as np
from numpy import * 
from copy import deepcopy
import matplotlib.pyplot as plt
from numpy import linspace
from math import pi, sin, factorial

def polyInter_coeff(fs,xs):
    fsc = deepcopy(fs)
    xsc = deepcopy(xs)
    xsc = transpose(xsc)
    length = len(xsc)
    A = matrix(zeros((length,length)))+ xsc
    for i in range(length):
        for j in range(length):
            A[i,j]= A[i,j]**(length-j-1)
    lists= fsc.tolist()[0]
    a_s = linalg.solve(A,lists)
    #print 'as = ',a_s
    return a_s
def Polynom(a_s,x):
    poly = 0.0
    count = len(a_s)-1
    for ii in range(len(a_s)):
        poly = poly + a_s[ii]*pow(x,count)
        count -= 1
    return poly


def lagrange(xs,x,i):
    lagr = 1.0
    for j in range(len(xs)):
        if i != j:
            lagr *= (x-xs[j])/(xs[i]-xs[j])
            #print 'lag ', lagr
    return lagr

def interLagrange(xs,fs,x):
        poly = 0.0
        for i in range(len(xs)):
            lagr = lagrange(xs, x, i)
            poly += fs[i]*lagr
        return poly

def div_diff(xs, starting_values):
    """Calculate the Newton dividing differences.

    Return the "top layer" which are the f(x0,...,xi) values
    """

    n = len(xs)
    fprev = starting_values
    fhist = []
    ret = []

    for k in range(1, n+1):
        fhist.append(fprev)
        ret.append(fprev[0])
        fprev = [(fprev[i+1] - fprev[i]) / (xs[k+i] - xs[i]) for i in range(n-k)]

    return (ret, fhist)

def horner(x, xs, diffs):
    n = len(xs)
    s = diffs[n-1]
    for k in range(n-2, -1, -1):
        s = s * (x - xs[k]) + diffs[k]
        
    return s

def newton_inter(xs, ys, x):
    """Interpolate f such that ys = f(xs) at the point x. """

    if len(xs) != len(ys):
        raise ValueError, "len(xs) != len(ys)"

    # Let's start with the dividing differences
    (diffs, _) = div_diff(xs, ys)
    
    # and use the Horner schema to interpolate
    return horner(x, xs, diffs)

class CubicSpline:
    def __init__(self, xs, ys):
        if len(xs) != len(ys):
            print 'xs ', xs
            print 'fs',ys
            raise ValueError, "len(xs) != len(ys)"

        self.xs = xs
        self.ys = ys
        self.moments()

    def moments(self):
        xs = self.xs
        ys = self.ys

        n = len(xs) - 1 # degrees of freedom

        mu = [0] * (n+1)
        for i in range(1, n):
            mu[i] = (xs[i] - xs[i-1])/(xs[i+1] - xs[i-1])
            mu[n] = (xs[n] - xs[n-1])/(xs[n] - xs[1] - xs[0] - xs[n-1])

        A = np.matrix(np.zeros((n-1, n-1)))
        #print "n", n, "A", A
        for i in range(n-1):
            print "i", i
            if i != 0:
                A[i,i-1] = mu[i+1]
            A[i,i] = 2
            if i != n-2:
                A[i,i+1] = 1.0 - mu[i+1]

        b = np.matrix(np.zeros((n-1, 1)))
        (d, dd) = div_diff(xs, ys)
        print "d", d, "dd", dd
        for i in range(0, n-1):
            b[i,0] = 6.0*dd[2][i]

        [self.M] = np.linalg.solve(A, b).flatten().tolist()

    def __call__(self, x):
        xs = self.xs
        ys = self.ys
        n = len(xs)
        # Add x_{n+1} and y)_{n+1} to get the edge cases
        #xs2 = xs + [xs[-1] + xs[1] + xs[0]]
        #ys2 = ys + [ys[-1]]

        M = [0.0] + self.M + [0.0]
        #M.insert(0, 0.0)
        #M.append(0.0)

        h = [xs[i] - xs[i-1] for i in range(n)]

        C = [0] * (n)
        for i in range(n):
            first = (ys[i] - ys[i-1]) / h[i]
            second = (h[i]/6.0) * (M[i] - M[i-1])
            C[i] = first - second

        # def C(i):
            # first = (ys[i] - ys[i-1]) / h[i]
            # second = (h[i]/6.0) * (M[i] - M[i-1])
            # return first - second

        D = [0] * n
        for i in range(n):
            first = (ys[i] + ys[i-1]) / 2.0
            second = ((h[i]**2)/12.0) * (M[i] + M[i-1])
            D[i] = first - second

        # def D(i):
        #     first = (ys[i] + ys[i-1]) / 2.0
        #     second = ((h[i]**2)/12.0) * (M[i] + M[i-1])
        #     return first - second

        def s(x, i):
            first = M[i-1] * (((xs[i] - x)**3) / (6*h[i]))
            second = M[i] * (((x - xs[i-1])**3) / (6*h[i]))
            third = C[i] * (x - ((xs[i-1]+xs[i])/2.0))
            fourth = D[i]
        #print 'first', first, 'second', second, 'third', third, 'fourth', fourth
            return first + second + third + fourth

        # find what interval x is in
        i = -1
        for ii in range(1, n):
            if x >= xs[ii-1] and x <= xs[ii]:
                i = ii
                break

        if i == -1:
            raise ValueError, ("No such interval for %f" % x)

        return s(x, i)




def main():
    #Eingabe
    xs = matrix(range(1,14),dtype= float)
    fs =matrix ('3.4 4.1 2.2 2.9 2.9 3.1 3.6 4.1 3.9 1.1 0.7 1.1 2.1', dtype= float)
    
    def polyInter(anzahl):
        print '### 5.1.1 ###'
        anzahl = 14/anzahl
        test_xs = np.arange(0,14,anzahl)
        a_s = polyInter_coeff(fs,xs)
        lpolyInter = []
        for jj in test_xs:
            poly = Polynom(a_s,jj)
            lpolyInter.append(poly)
        print 'Egebnis von Polynominterpolation: ', lpolyInter
        #plt.plot(test_xs,lpolyInter,'.-r')
        return (test_xs,lpolyInter)
        
    def polyLag(anzahl):
        print '### 5.1.2 ###'
        test_xs = np.arange(0,14,14/anzahl)
        list_lagr = []
        for jj in test_xs:
            p = interLagrange(xs.tolist()[0], fs.tolist()[0], jj)
            list_lagr.append(p)
        print 'Ergebnis von Lagrange-Interpolation', list_lagr
        #plt.plot(test_xs,list_lagr,'o-b')
        return (test_xs,list_lagr)
        
    def polyNewt(anzahl):
        print '### 5.1.3 ###'
        l = []
        anzahl = 14/anzahl
        test_xs = np.arange(0,14,anzahl)
        for k in test_xs:
            poly = newton_inter(xs.tolist()[0], fs.tolist()[0], k)
            l.append(poly)
        print 'Ergebnis von Newton-Interpolation' ,l
        #plt.plot(test_xs,l,'+-g')
        #plt.legend((u'Gegeben',))
        #plt.title(ur'Interpolation von $\sin(2\pi x)$')
        return (test_xs,l)

    def poly_cspline(anzahl):
        print '### 5.1.4 ###'
        #Es = []
        lists = []
        x_list = []
        for n in range(2, 15,2):
            xss = xs.tolist()[0][0:n]
            fss = fs.tolist()[0][0:n]
            cs = CubicSpline(array(xss), array(fss))
            #print "M", M
            #errors = []
            points = linspace(xss[n-1], xss[n], anzahl, endpoint=False)
            x_list.extend(points)
            for p in points:
                lists.append(cs(p))
                #errors.append(abs(cubic_splines(M, xs, ys, p) - sin(2*pi*p)))
            #Es.append((n, max([e.max() for e in errors])))
            print 'Ergebnis von Cubic-Splines: ', lists
        return (x_list,lists)
    
    def compare_speed(anzahl):
        print '### 5.1.5 ###'
        time_polyInt = []
        time_lagr = []
        time_newton = []
        time_spline = []
        for anz in np.arange(anzahl):
            
            start = time.clock()
            polyInter(anz)
            elapsed = (time.clock() - start)
            time_polyInt.append(elapsed)
            
            start = time.clock()
            polyLag(anz)
            elapsed = (time.clock() - start)
            time_lagr.append(elapsed)
            
            start = time.clock()
            polyNewt(anz)
            elapsed = (time.clock() - start)
            time_newton.append(elapsed)
            
            start = time.clock()
            poly_cspline(anz)
            elapsed = (time.clock() - start)
            time_spline.append(elapsed)
            
        plt.plot(np.arange(anzahl),time_polyInt,'.r')
        plt.plot(np.arange(anzahl),time_lagr,'.g')
        plt.plot(np.arange(anzahl),time_newton,'.b')   
        plt.legend((u'Polynominterpolation',u'Lagrange-Interpolation',u'Newton-Interpolation'))
    
    
    plt.plot(xs.tolist()[0],fs.tolist()[0])
    '''
    polyInter(30.)
    polyLag(30.)
    polyNewt(30.)
    # hier ist die Anzahl der x-werte zwischen x{i-1} und x{i}
    poly_cspline(20)
    plt.xlabel(u'Interpolationspunkte')
    plt.ylabel(u'P(x)')
    '''
    # hier ist die Anzahl der x-werte zwischen x{i-1} und x{i}
    
    plt.plot(zip(*poly_cspline(20)))
    #compare_speed(300.)
    '''
    Aufgabe 5.1.5:
    Vergleich der Zeitmessung zwischen vier Verfahren bei der Anzahl der Elementen n = 0..500 (Siehe Diagramm):
    
    1) Normale Polynominterpolation: am schnellstens, sehr leichter Anstieg der Zeit
    2) Lagrange-Interpolation: am langsamstens, schneller Anstieg der Zeit als die anderen
    3) Newton-Interpolation: zweit schnellstes Verfahren nach Poly.interp.
    4) Kubisches Spline: TODO
    '''
    
if __name__ == '__main__':
    main()
    plt.show()
