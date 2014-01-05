#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
import math

def diff_Vx(t,vx,vy):
    return -0.004*np.sqrt(vx**2 + vy**2)*vx
def diff_Vy(t,vx,vy):
    return -0.004*np.sqrt(vx**2 + vy**2)*vy-9.8

def rungeKutta_9_1(Fx,Fy,a,b,h,vx_0,vy_0,x0,y0):
    vXn = vx_0
    vYn = vy_0
    xn = x0
    yn = y0
    tn = a
    n = int(round((b-a) / h))
    #print("[RUNGE KUTTA] tn = %f yn = %s"%(tn, yn))
    
    ts = [a]
    vx = [vx_0]
    vy = [vy_0]
    x = [x0]
    y = [y0]
    def F_x(t,v):
        return v;
    def F_y(t,v):
        return v;
    
    while n > 0:
        #(t,k,vx,vy)
        
        k1 = F_y(tn,vYn)
        k2 = F_y(tn + 0.5*h ,vYn + h*0.5*k1)
        k3 = F_y(tn + 0.5*h ,vYn + h*0.5*k2)
        k4 = F_y(tn + h , vYn + h*k3)
        
        yn = yn + h/6.0*(k1 + 2*k2 + 2*k3 + k4)
        
        k1_ = F_x(tn,vXn)
        k2_ = F_x(tn + 0.5*h ,vXn + h*0.5*k1_)
        k3_ = F_x(tn + 0.5*h ,vXn + h*0.5*k2_)
        k4_ = F_x(tn + h , vXn + h*k3_)
        
        xn = xn + h/6.0*(k1_ + 2*k2_ + 2*k3_ + k4_)
        
        k1x = Fx(tn,vXn,vYn)
        k2x = Fx(tn + 0.5*h , vXn + h*0.5*k1x,vYn)
        k3x = Fx(tn + 0.5*h , vXn + h*0.5*k2x,vYn)
        k4x = Fx(tn + h , vXn + h*k3x,vYn)
        
        k1y = Fy(tn,vXn,vYn)
        k2y = Fy(tn + 0.5*h , vXn,vYn + h*0.5*k1y)
        k3y = Fy(tn + 0.5*h , vXn,vYn + h*0.5*k2y)
        k4y = Fy(tn + h , vXn, vYn + h*k3y)
        
        vXn = vXn + h/6.0*(k1x + 2*k2x + 2*k3x + k4x)
        vYn = vYn + h/6.0*(k1y + 2*k2y + 2*k3y + k4y)
        
        '''
        xn = xn + vXn*h
        yn = yn + vYn*h
        vXn = vXn + h*Fx(tn,vXn,vYn)
        vYn = vYn + h*Fy(tn,vXn,vYn)
        '''
        if xn >= 200.0 and xn < 200.1 and yn > 0.0 and (yn+tn) < 51.0 and (yn+tn) >= 50.0:
            print("tn = %f vxn = %f vyn = %f xn = %f yn = %f"%(tn,vXn,vYn,xn,yn))
        tn += h
        
        n -= 1
        
        # Ergebnisse speichern
        ts.append(tn)
        vx.append(vXn)
        vy.append(vYn)
        x.append(xn)
        y.append(yn)
        #print("[RUNGE KUTTA] tn = %f yn = %s"%(tn, yn))
        
    return ts,vx,vy,x,y

def run_9_1_1():
    def alpha_func(alpha):
        print 'Alpha = ', alpha
        vx_0 = abs(100 * math.cos(alpha)) # satz der pythagoras
        vy_0 = abs(100 * math.sin(alpha))# satz der pythagoras
        x0 = 0.0
        y0 = 0.0
        a = 0.0
        b = 15
        h = 0.001 
        ts,vx,vy,x,y = rungeKutta_9_1(diff_Vx,diff_Vy,a,b,h,vx_0,vy_0,x0,y0)
        plt.plot(x,y)
    for i in np.arange (1,30,0.1):
        #alpha_func(5)
        alpha_func(i)
        #alpha_func(50)
        #alpha_func(70)
    plt.show()
        
#Ergebnis fuer 9.1.2 ohne Annäherung:
#Alpha =  13
#tn = 3.463000 vxn = 39.456931 vyn = -6.248466 xn = 200.038937 yn = 46.859919
#Alpha =  24.7
#tn = 3.463000 vxn = 39.456931 vyn = -6.248466 xn = 200.038937 yn = 46.859919

    
def nst_f(alpha):
    #Alpha = 26
    #tn = 3.462000 vxn = 39.474342 vyn = -6.275243 xn = 200.074154 yn = 46.686820
    return -6.275243*3.462*math.sin(alpha) - 1/2*9.8*3.462**2 - 50.0 + -6.275243*3.462

def fallschirm_hoehe(t):
    #vp = 1
    return 50.0 - t
        
def next_point(f, x):
    """ Return the next point for Newton """
    deriv = derivative(f, x)
    return x - (f(x)/deriv)

def newton(f, x0, delta):
    """
    Return the zero with `delta` precision
    """

    last_point = x0
    point = next_point(f, x0)

    n = 0
    while not abs(last_point - point) < delta:
        last_point = point
        point = next_point(f, last_point)
        print point
        n += 1

    return point, n

def sekante(f,x0,delta):
    last_point = 0.0
    point = x0

    n = 0
    while not abs(last_point - point) < delta:
        tmp_point = point - ((point - last_point)/(f(point)-f(last_point))) * f(point)
        last_pont = point
        point = tmp_point
        print point
        n += 1

    return point, n

if __name__ == '__main__':
    #point,n = newton(nst_f,5.0,0.00001)
    run_9_1_1()
    #print fallschirm_hoehe(2.527100)
