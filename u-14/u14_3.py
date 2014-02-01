# -*- coding: utf-8 -*-
'''
@author: Carlos Matin Nieto, Tu Tran
'''
import numpy as np
import matplotlib.pyplot as plt

def rw(N,d=1):
    res = np.zeros((N,d))
    for i in range(N):
        r = np.random.random_integers(0,1)
        if r == 0:
            r = -1
        res[i][np.random.random_integers(d)-1] = r
    return res

def run():
    def predict_return(N,d,Pfade):
        anz_zurueck = 0
        for _ in range(Pfade):
            walks = rw(N,d)
            weite = np.zeros(d)
            for i in range(N):
                weite = weite + walks[i]
                if np.sum(np.abs(weite)) == 0:
                    anz_zurueck += 1
                    break
        print np.float(anz_zurueck)/np.float(Pfade) 
        return anz_zurueck/np.float(Pfade) 
    print '1D'
    d1 = [predict_return(np.int(np.exp(i)),1,200) for i in range(15)]
    #print d1
    print '2D'
    d2 = [predict_return(np.int(np.exp(i)),2,200) for i in range(15)]
    print '3D'
    d3 = [predict_return(np.int(np.exp(i)),3,200) for i in range(15)]
    plt.plot(range(15),d1,label = 'Prediction in 1D')
    plt.plot(range(15),d2,label = 'Prediction in 2D')
    plt.plot(range(15),d3,label = 'Prediction in 3D')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    run()
