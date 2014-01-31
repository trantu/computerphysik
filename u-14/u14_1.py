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

def predictedX(N,d=1):
    x0 = 0
    walks = rw(N,d)
    weite = np.zeros(d)
    var = np.zeros((N,d))
    for i in range(N):
        weite += walks[i]
        var[i] = (weite-x0)**2
    return np.sum(var) / N
def analytisch_X(N):
    return N
def u14_1_1():
    Nmax = 250
    laeufe = [200, 500, 1000]
    for i  in laeufe:
        Xs = []
        Xs_analytisch = []
        Xs_standard = []
        Xs_st_analy = []
        for j in range(1,Nmax):
            tmp = 0.0
            for _ in range(i):
                tmp += predictedX(j)
            
            Xs_analytisch.append(j)
            Xs.append(tmp/i)
            Xs_standard.append(np.sqrt(tmp/i))
            Xs_st_analy.append(np.sqrt(j))
        plt.plot(range(1,Nmax),Xs_analytisch,label = 'Xs_analytisch')
        plt.plot(range(1,Nmax),Xs,label = 'Xs')
        plt.plot(range(1,Nmax),Xs_standard,label = 'Xs_standard')
        plt.plot(range(1,Nmax),Xs_st_analy,label = 'Xs_st_analy')
        plt.title('Laeufe = '+str(i))
        plt.legend()
        plt.show()
    
    D = (sum([predictedX(Nmax) for _ in range(500)])/500) / (2.0*Nmax) # Nach der Gleichung in der Aufgabenstellung
    print D
    '''
    Diffusionskonstante = 0.251671552
    '''
def u14_1_2():
    D_2dim = (sum([predictedX(250,2) for _ in range(200)])/200) / (2.0*250)
    print 'Diffusionskons in 2D', D_2dim
    D_3dim = (sum([predictedX(250,3) for _ in range(200)])/200) / (2.0*250)
    print 'Diffusionskons in 3D', D_3dim
    '''
    Diffusionskons in 2D: 0.26368264
    Diffusionskons in 3D: 0.24547456
    also bleibt beinah unverändert
    '''
    laeufe = 2000
    Nmax = 500
    predict_2dim = sum([predictedX(Nmax, 2) for _ in range(laeufe)])/laeufe
    print 'Erwartungswert für 2-dimensional',predict_2dim
    predict_3dim = sum([predictedX(Nmax, 3) for _ in range(laeufe)])/laeufe
    print 'Erwartungswert für 3-dimensional',predict_3dim
    '''
    Erwartungswert für 2-dimensional 240.444554
    Erwartungswert für 3-dimensional 253.699132
    '''
    def entfernungen(N,d=1):
        walks = rw(N,d)
        weite = np.zeros(d)
        var = []
        for i in range(N):
            weite += walks[i]
            var.append(np.sum(weite))
        return var
    #for _ in range(200):
    plt.hist(entfernungen(250),normed=True, color="g",label = '1D')
    plt.hist(entfernungen(250,2),normed=True, color="y",label = '2D')
    plt.hist(entfernungen(250,3),normed=True, color="b",label = '3D')
    plt.legend()
    plt.show()
if __name__ == '__main__':
    u14_1_1()
    u14_1_2()
