'''
Created on Nov 14, 2013

@author: calosnm,tran
'''
# encoding: utf-8

import math
import copy
from numpy import * 
from copy import deepcopy

def tausch_Zeilen(m, a, b):
    temp = matrix(m[a])
    m[a] = m[b]
    m[b] = temp
    
def mround ( x,N ) :
    if ( x== 0.0 ) :
        return x
    return round ( x, int ( N - math.ceil ( math.log10 ( abs ( x )))))

def eliminable(A,b):
    zeilen,spalten = shape(A)
    if (zeilen != spalten or zeilen == 0 or spalten == 0):
        raise ValueError("ungleiche Anzahl zwischen Zeilen und Spalten")
    for i in range(zeilen):
        if mean(A[i,:]) == 0.0:
            raise ValueError("eine Zeile Null")
        
def konvergiert(xk,xk_1):
    diff = xk_1-xk
    skalar = 0.0
    zeilen,_ = shape(xk)
    genauigkeit = 0.5*10**(-4)
    for i in range(zeilen):
        skalar = skalar + diff[i]**2
    if(sqrt(skalar)>genauigkeit):
        return 0
    else:
        return 1

def gauss_seidel_alg(A,b):
    F = deepcopy(A)
    d = deepcopy(b)
    L = deepcopy(A)
    U = deepcopy(A)
    D = deepcopy(A)
    zeilen,spalten = shape(A)
    #L, U, D bekommen
    for i in range(spalten):
        for j in range(zeilen):
            if i != j:
                D[i,j] = 0
                if i < j:
                    L[i,j] = 0
                else:
                    U[i,j] = 0
            else:
                U[i,j] = 0
                L[i,j] = 0
                
    #gauss Seidel Teil
    L_D = D+L
    L_D_inverse = L_D**(-1)
    x_i = matrix('0;0;0')
    x_i_plus = matrix('0;0;0')
    nb_iteration = 0
    while 1:
        x_i_plus = L_D_inverse*(b-U*x_i)
        if konvergiert(x_i, x_i_plus):
            
            break
        x_i = x_i_plus
        nb_iteration = nb_iteration +1
    
    for k in range(zeilen):
        x_i[k] = mround(x_i[k],4)
    print 'Loesung'
    print x_i
    print 'Anzahl der Iteration'
    print nb_iteration
    return x_i, nb_iteration
    
'''
A = matrix('16 3;7 -11',dtype = float)
b = matrix('11; 13',dtype = float)
c = matrix('1;1;0',dtype = float)
'''
A = matrix('3.5 3 0; -1 4 4; 0 3 4.5',dtype = float)
A[2,0] = 1/3.
A[0,2] = -0.5
b = matrix('7.5; -6.5; 1',dtype = float)
print 'Aufgabe 4.1.1'
gauss_seidel_alg(A, b)
#konvergiert(b, c)

#TODO Testvektor xk wird nach Abschluss der Iteration aktulisiert
def Jacobi(A,b):
    d = deepcopy(b)
    L = deepcopy(A)
    U = deepcopy(A)
    D = deepcopy(A)
    zeilen,spalten = shape(A)
    #L, U, D bekommen
    for i in range(spalten):
        for j in range(zeilen):
            if i != j:
                D[i,j] = 0
                if i < j:
                    L[i,j] = 0
                else:
                    U[i,j] = 0
            else:
                U[i,j] = 0
                L[i,j] = 0
                
    #jacobi teil
    D_inverse = D**(-1)
    x_i = matrix('0;0;0')
    x_i_plus = matrix('0;0;0')
    nb_iteration = 0
    while 1:
        x_i_plus = D_inverse*((-L-U)*x_i+d)
        if konvergiert(x_i, x_i_plus):
            break
        x_i = x_i_plus
        nb_iteration = nb_iteration +1
    
    for k in range(zeilen):
        x_i[k] = mround(x_i[k],4)
    print 'Loesung'
    print x_i
    print 'Anzahl der Iteration'
    print nb_iteration
    return x_i, nb_iteration

print 'Aufgabe 4.1.2'
Jacobi(A, b)
#Man sieht bei Jacobi gibt es mehr Iterationen

# Aufgabe 4.1.3: TODO
# Aufgabe 4.1.4:

