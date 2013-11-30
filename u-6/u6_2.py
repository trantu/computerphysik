'''
Created on Nov 29, 2013

@author: Carlos Matin Nieto, Tu Tran
'''
# -*- coding: utf-8 -*-
import time
import math
import copy
import numpy as np
from numpy import * 
from copy import deepcopy
import matplotlib.pyplot as plt
from numpy import linspace
from math import pi, sin,cos, factorial,exp

def f(a0,a1,a2,a3,t):
    return (exp(-a0 * t )) * (a2*sin(a1*t) + a3*cos(a1*t))
def g(a0,a1,a2,t):
    return (exp(-a0 * t)) * a2 * sin(a1 * t)

#def gauss_newton(jacob):
def jacob_f(xs,a0,a1,a2,a3):
    
    def f_a0(x):
        return (-x*exp(-a0*x))* (a2*sin(a1*x) + a3*cos(a1*x))         
    def f_a1(x):
        return (x*exp(-a0*x))* (a2*cos(a1*x) - a3 * sin(a1*x))
    def f_a2(x):
        return (exp(-a0*x))* sin(a1*x)
    def f_a3(x):
        return (exp(-a0*x))* cos(a1*x)
    jmatrix_f = np.matrix(np.zeros((len(xs), 4)))
    for i in range(len(xs)):
        jmatrix_f[i] = [f_a0(xs[i]), f_a1(xs[i]), f_a2(xs[i]), f_a3(xs[i])]
        
    return jmatrix_f
    
def jacob_g(xs,a0,a1,a2):
    
    def g_a0(x):
        return (-x*exp(-a0*x))* a2*sin(a1*x)         
    def g_a1(x):
        return (x*exp(-a0*x))* a2*cos(a1*x)
    def g_a2(x):
        return (exp(-a0*x))*sin(a1*x)
    
    jmatrix_g = np.matrix(np.zeros((len(xs), 4)))
    for i in range(len(xs)):
        jmatrix_g[i] = [g_a0(xs[i]), g_a1(xs[i]), g_a2(xs[i])]
    return jmatrix_g
