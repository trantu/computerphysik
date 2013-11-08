#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

def tausch_Zeilen(m,z1,z2):
    for i in range(4):
        wert1 = m[z1,i]
        m[z1,i] = m[z2,i]
        m[z2,i] = wert1
    return m
def gauss_alg(F, d):
    zeilen,spalten = shape(F)
    # nicht moeglich zu eliminieren
    if (zeilen != spalten | zeilen == 0 | spalten == 0):
        raise ValueError("ungleiche Anzahl zwischen Zeilen und Spalten")
    # axis = 1 -> vertical
    F = concatenate((F, d.T), axis=1)
    print F
    pivot = F[0,0]
    #Zeilenvertauschen
    if(pivot == 0 and (F[1,0] != 0 or F[2,0] !=0)):
        if(F[1,0] != 0 ):
            F = tausch_Zeilen(F, 0, 1)
        else:
            F = tausch_Zeilen(F, 0, 2)
    pivot = F[0,0]    
    if(F[1,0] != 0):
        alp = alpha(F[1,0], pivot)
        for i in range(4):
            F[1,i]= (alp* F[0,i])+F[1,i]
    if(F[2,0] != 0):
        alp = alpha(F[2,0], pivot)
        for i in range(4):
            F[2,i]= (alp* F[0,i])+F[2,i]
    # Addition mit 2. Zeile
    pivot = F[1,1]
    if(pivot == 0 and F[2,1] != 0):
            F = tausch_Zeilen(F, 1, 2)
            print 'tausch'
    pivot = F[1,1]
    if(F[2,1] != 0):
        alp = alpha(F[2,1], pivot)
        print alp
        for i in range(1,4):
            F[2,i]= (alp* F[1,i])+F[2,i]
    return F
'''
Test
'''
m = matrix ('1 -2 -1; 36000 2 0; -2 1400 1',dtype = float);
d = matrix ('3 72002 1399',dtype = float);
'''
Ausgabe:
[[  1.00000000e+00  -2.00000000e+00  -1.00000000e+00   3.00000000e+00]
 [  0.00000000e+00   7.20020000e+04   3.60000000e+04  -3.59980000e+04]
 [  0.00000000e+00   0.00000000e+00  -6.98980612e+02   2.10294183e+03]]
'''
print gauss_alg(m, d)


