'''
Created on Nov 14, 2013

@author: calosnm,tran
'''
# encoding: utf-8

import math
import copy
from numpy import * 
from copy import deepcopy
import matplotlib.pyplot as plt

def tausch_Zeilen(m, a, b):
    temp = matrix(m[a])
    m[a] = m[b]
    m[b] = temp
    
def mround ( x,N ) :
    if ( x== 0.0 ) :
        return x
    return round ( x, int ( N - math.ceil ( math.log10 ( abs ( x )))))
        
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
def solve (F, d):
    A = deepcopy(F)
    b = deepcopy(d)
    zeilen,_ = shape(A)
    for i in range(zeilen):
        for j in range(zeilen):
            if j == i:
                A[i,j] = b[i]/A[i,j]
                for ii in range(i+1,zeilen):
                    A[ii,j] = A[ii,j]*A[i,j]
                break
            b[i] = b[i] - A[i,j]
            A[i,j] = 0.0
    #print 'Result: [X1 X2 X3]'
    #print diag(A)
    return transpose(matrix(diag(A)))

def gauss_seidel_alg(A,b,x_ini):
    print '### Gauss-Seidel ###'
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
    #L_D_inverse = L_D**(-1)
    x_i = deepcopy(x_ini)
    x_i_plus = deepcopy(x_ini)
    nb_iteration = 0
    while 1:
        nb_iteration = nb_iteration +1
        #vorwaerst einsetzen
        x_i_plus = solve(L_D,-U*x_i+d)
        if konvergiert(x_i, x_i_plus) or nb_iteration > 500:            
            break
        x_i = x_i_plus
    for k in range(zeilen):
        x_i[k] = mround(x_i[k],4)
    print 'Loesung: '
    print x_i
    print 'Anzahl der Iteration: '
    print nb_iteration
    return x_i, nb_iteration
    
'''
x_ini = matrix('0;0')
A = matrix('16 3;7 -11',dtype = float)
b = matrix('11; 13',dtype = float)
'''
A = matrix('3.5 3 0; -1 4 4; 0 3 4.5',dtype = float)
A[2,0] = 1/3.
A[0,2] = -0.5
b = matrix('7.5; -6.5; 1',dtype = float)

x_ini = matrix('0;0;0')
print 'Aufgabe 4.1.1'
gauss_seidel_alg(A, b,x_ini)
#solve(A,b)
#konvergiert(b, c)


#TODO Testvektor xk wird nach Abschluss der Iteration aktulisiert
def Jacobi(A,b,x_ini):
    print '#### Jacobi ####'
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
    #D_inverse = D**(-1)
    x_i = deepcopy(x_ini)
    x_i_plus = deepcopy(x_ini)
    nb_iteration = 0
    while 1:
        nb_iteration = nb_iteration +1
        x_i_plus = solve(D,((-L-U)*x_i+d))
        if konvergiert(x_i, x_i_plus) or nb_iteration > 500:
            break
        x_i = x_i_plus
    
    for k in range(zeilen):
        x_i[k] = mround(x_i[k],4)
    print 'Loesung:'
    print x_i
    print 'Anzahl der Iteration:'
    print nb_iteration
    return x_i, nb_iteration

print 'Aufgabe 4.1.2:'
Jacobi(A, b,x_ini)
#Man sieht in diesem Beispiel bei Jacobi gibt es mehr Iterationen im Vergleich zu Gauss-Seidel

# Aufgabe 4.1.3: Die Konvergenz ist nur garantiert, wenn die Matrix entweder diagonaldominant, symmetrisch oder positiv definit ist.

print 'Aufgabe 4.1.4:'
A4 = matrix('5 3 -1 2; -3 7 6 -2; 4 4 3 -3; -5 2 2 4',dtype = float)
b4 = matrix('8; 1; 7; 2',dtype = float)
x_ini4 = matrix('0;0;0;0')
gauss_seidel_alg(A4, b4, x_ini4)
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
    list_iter = []
    for w in relax:
        #L_D = (1/w)*D+L
        #L_D_inverse = L_D**(-1)
        x_i = deepcopy(x_ini)
        x_i_plus = deepcopy(x_ini)
        nb_iteration = 0
        while 1:
            nb_iteration = nb_iteration +1
            x_i_plus = solve(L+(1/w)*D,((1/w)*D-D-U)*x_i+b)
            if konvergiert(x_i, x_i_plus) or nb_iteration > 250:               
                break
            x_i = x_i_plus
        
        for k in range(zeilen):
            x_i[k] = mround(x_i[k],4)
        list_iter.append(nb_iteration)
        print 'Loesung: '
        print x_i
        print 'Anzahl der Iteration: '
        print nb_iteration
        print 'Relax w: %f' %w
    return x_i, nb_iteration,list_iter,relax

nr_of_relax = 50
print 'Aufgabe 4.1.5:'
x_i, nb_iteration,list_iter,relax = gauss_seidel_alg_relax(A4, b4, x_ini4, nr_of_relax)
# bei w = 0.88 ist werden 39 Iterationen benoetigt, also am wernigstens. Ab w > 1 scheinen die Werte nicht konvergieren

plt.plot(relax,list_iter)
plt.xlabel("$\omega$ - Werte")
plt.ylabel("Anzahl der Iteration")
plt.show()
