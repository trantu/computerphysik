'''
Created on Nov 8, 2013

@author: tran, carlosmn
'''
# encoding: utf-8
import math
import copy
from numpy import * 
from copy import deepcopy


def mround ( x,N ) :
    if ( x== 0.0 ) :
        return x
    return round ( x, int ( N - math.ceil ( math.log10 ( abs ( x )))))
def alpha(a,pivot):
    return -1*(a/pivot)

def tausch_Zeilen(m, a, b):
    temp = matrix(m[a])
    m[a] = m[b]
    m[b] = temp
def eliminable(A,b):
    zeilen,spalten = shape(A)
    if (zeilen != spalten or zeilen == 0 or spalten == 0):
        raise ValueError("ungleiche Anzahl zwischen Zeilen und Spalten")
    for i in range(zeilen):
        if mean(A[i,:]) == 0.0:
            raise ValueError("eine Zeile Null")
    return 1
#rueckwaerts loesen
def solve (F, d,N):
    A = deepcopy(F)
    b = deepcopy(d)
    zeilen,_ = shape(A)
    rng = range(zeilen)
    rng.reverse()
    for i in rng:
        for j in rng:
            if j == i:
                A[i,j] = mround(b[i]/A[i,j],N)
                for ii in range(i):
                    A[ii,j] = mround(A[ii,j]*A[i,j],N)
                break
            b[i] = mround(b[i] - A[i,j],N)
            A[i,j] = 0.0
    print 'Result: [X1 X2 X3]'
    print diag(A)
    return diag(A)

def gauss_alg_mround(A, b, N):
    F = deepcopy(A)
    d = deepcopy(b)
    zeilen,spalten = shape(F)
    d = d.T
    for s in range(spalten):
        #Ueberpruefen, ob eliminierbar
        eliminable(F, d)
        pivot = F[s,s]
        #Zeilenvertauschen
        if(pivot == 0):
            for ii in range(s+1,zeilen):
                if(F[ii,s] != 0):
                    tausch_Zeilen(F, s, ii)
                    tausch_Zeilen(d, s, ii)
                    break
        pivot = F[s,s]
        for zz in range(s+1,zeilen):  
            #addition der Zeilen  
            if F[zz,s] != 0:
                alp = alpha(F[zz,s], pivot)
                d[zz] = mround((alp* d[s])+d[zz], N)
                for i in range(s,spalten):
                    F[zz,i] = mround((alp * F[s,i]) + F[zz,i], N)
    print 'Matrix A'
    print F
    print 'Rechte Seite b'
    print d
    return F,d
'''
# 3.1.2: Im Vergleich zu dem Verfahren ohne mround mit N = 4 sieht man mit der Rundenfunktion, dass die Zahlen an der 4. Stelle nach dem Komma
gerundet werden 
'''
# Suche max Element in einer Spalte und gebe die Zeilennummer zuruek
def max_pivot(spalte,zeilennummer):
    a,_ = shape(spalte)
    zeile = 0
    maxi = 0.0
    for i in range(zeilennummer,a):
        if abs(spalte[i]) > abs(maxi):
            maxi = spalte[i]
            zeile = i
    return zeile

def gauss_alg_mround_best_pivot(A, b, N):
    F = deepcopy(A)
    d = deepcopy(b)
    zeilen,spalten = shape(F)
    d = d.T
    for s in range(spalten):
        eliminable(F, d)
        getauscht = False
        if ~getauscht:
            zeile = max_pivot(F[:,s],s)
            # tausche max pivot
            if zeile != s:
                tausch_Zeilen(F, s, zeile)
                tausch_Zeilen(d, s, zeile)
            getauscht = True
        #print F    
        pivot = F[s,s]
        for zz in range(s+1,zeilen):  
            #addition der Zeilen  
            if F[zz,s] != 0:
                alp = alpha(F[zz,s], pivot)
                d[zz] = mround((alp* d[s])+d[zz], N)
                for i in range(s,spalten):
                    F[zz,i] = mround((alp * F[s,i]) + F[zz,i], N)
    print 'Matrix A'
    print F
    print 'Rechte Seite b'
    print d
    return F, d
        
m = matrix ('1 -2 -1; 36000 2 0; -2 1400 1',dtype = float);
d = matrix ('3 72002 1399',dtype = float);

print '##### 3.1.2 #####'
# N = 20 ist nur symbolisch. Es soll hier nicht gerundet werden.
m1,d1 = gauss_alg_mround(m, d, 20)
solve(m1,d1,10)

print '##### 3.1.2 #####'
print 'matrix mit mround , N = 4'
m2,d2 = gauss_alg_mround(m, d, 4)
solve(m2,d2,4)
print '##### 3.1.3 #####'
print 'matrix mit mround , N = 4'
m3,d3 = gauss_alg_mround_best_pivot(m, d, 4)
solve(m3,d3,4)
