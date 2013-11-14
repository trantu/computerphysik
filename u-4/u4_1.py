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

def gauss_seidel_alg(A,b,x_ini):
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
    x_i = deepcopy(x_ini)
    x_i_plus = deepcopy(x_ini)
    nb_iteration = 0
    while 1:
        x_i_plus = L_D_inverse*(b-U*x_i)
        if konvergiert(x_i, x_i_plus) or nb_iteration > 500:
            
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
x_ini = matrix('0;0;0')
print 'Aufgabe 4.1.1'
gauss_seidel_alg(A, b,x_ini)
#konvergiert(b, c)

#TODO Testvektor xk wird nach Abschluss der Iteration aktulisiert
def Jacobi(A,b,x_ini):
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
    x_i = deepcopy(x_ini)
    x_i_plus = deepcopy(x_ini)
    nb_iteration = 0
    while 1:
        x_i_plus = D_inverse*((-L-U)*x_i+d)
        if konvergiert(x_i, x_i_plus) or nb_iteration > 500:
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
Jacobi(A, b,x_ini)
#Man sieht bei Jacobi gibt es mehr Iterationen

# Aufgabe 4.1.3: TODO
print 'Aufgabe 4.1.4:'
A4 = matrix('5 3 -1 2; -3 7 6 -2; 4 4 3 -3; -5 2 2 4',dtype = float)
b4 = matrix('8; 1; 7; 2',dtype = float)
x_ini4 = matrix('0;0;0;0')
gauss_seidel_alg(A4, b4, x_ini4)
print 'Jacobi:'
Jacobi(A4, b4, x_ini4)
# Die Gleichung laesst sich von Gauss-Seidel-Verfahren loesen, bei Jacobi kovergieren die Werte nicht trotz der vielen Iteration.
def gauss_seidel_alg_relax(A,b,x_ini,nr_of_relax):
    relax = [(2./nr_of_relax)*i for i in range(1,nr_of_relax)]
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
    x_i = deepcopy(x_ini)
    x_i_plus = deepcopy(x_ini)
    nb_iteration = 0
    for w in relax:
        L_D = (1/w)*D+L
        L_D_inverse = L_D**(-1)
        while 1:
            x_i_plus = L_D_inverse*(b-((1/w)*D-D-U)*x_i)
            if konvergiert(x_i, x_i_plus) or nb_iteration > 100000:               
                break
            x_i = x_i_plus
            nb_iteration = nb_iteration +1
        
        for k in range(zeilen):
            x_i[k] = mround(x_i[k],4)
        print 'Loesung'
        print x_i
        print 'Anzahl der Iteration'
        print nb_iteration
        print 'Relax w: %f' %w
    return x_i, nb_iteration

nr_of_relax = 30
print 'Aufgabe 4.1.4:'
gauss_seidel_alg_relax(A4, b4, x_ini4, nr_of_relax)
# w<=1 dann konvergieren die Werte
