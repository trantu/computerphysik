#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## Carlos Martín, Tran Tu

import numpy as np
import matplotlib.pyplot as plt

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
    
    
    while n > 0:
        #(t,k,vx,vy)
        xn = xn + vXn * h
        yn = yn + vYn * h
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
    vx_0 = 100.0
    vy_0 = 100.0
    x0 = 0.0
    y0 = 0.0
    a = 0.0
    b = 70
    h = 0.01
    ts,vx,vy,x,y = rungeKutta_9_1(diff_Vx,diff_Vy,a,b,h,vx_0,vy_0,x0,y0)
    plt.plot(ts,vx)
    plt.show()
    plt.plot(ts,vy)
    plt.show()
    plt.plot(x,y)
    plt.show()
if __name__ == '__main__':
    run_9_1_1()
