'''
Created on Nov 8, 2013

@author: tran, carlosmn
'''
# encoding: utf-8
import math
import copy
from numpy import * 
def mround ( x,N ) :
    if ( x==0 ) :
        return x
    return round ( x, int ( N - math.ceil ( math.log10 ( abs ( x )))))

def alpha(a,pivot):
    return -1*(a/pivot)

def swap(m, a, b):
    temp = matrix(m[a])
    m[a] = m[b]
    m[b] = temp

# first, d is ignored
def gauss_alg(F, d):
    zeilen,spalten = shape(F)
    # nicht moeglich zu eliminieren
    if zeilen != spalten or zeilen == 0 or spalten == 0:
        raise ValueError("ungleiche Anzahl zwischen Zeilen und Spalten")

    # axis = 1 -> vertical
    #F = concatenate((F, d.T), axis=1)
    for s in range(spalten):
        for z in range(s,zeilen):
            pivot = F[z,z]
            #Zeilenvertauschen
            if pivot == 0:
                for ii in range(z+1,zeilen):
                    if F[ii,s] != 0:
                        swap(F, z, ii)
                        break # back to second loop
                    
            pivot = F[z,z]
            #print pivot
            for zz in range(z+1,zeilen):    
                if F[zz,s] != 0:
                    #print F[z,:]
                    #print F[z+1,:]
                    alp = alpha(F[zz,s], pivot)
                    F[zz] = (alp * F[z]) + F[zz]
    return F
'''
Test
'''
#m = matrix('1 2 3; 1 1 1; 3 3 1')
#d = matrix('2 2 0')
m = matrix ('1 -2 -1; 36000 2 0; -2 1400 1',dtype = float);
d = matrix ('3 72002 1399',dtype = float);
'''
Ausgabe:
[[  1.00000000e+00  -2.00000000e+00  -1.00000000e+00   3.00000000e+00]
 [  0.00000000e+00   7.20020000e+04   3.60000000e+04  -3.59980000e+04]
 [  0.00000000e+00   0.00000000e+00  -6.98980612e+02   2.10294183e+03]]
'''
print gauss_alg(m, d)


