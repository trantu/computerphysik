# -*- coding: utf-8 -*-
'''
@author: Carlos Matin Nieto, Tu Tran
'''

import numpy as np
import matplotlib.pyplot as plt

def spielen(M,time):
    entwicklung = []
    while M > 0 and time > 0:
        entwicklung.append(M)
        if np.random.random_integers(0,1):
            M = M + 1
        else: M = M - 1 
        time -= 1
    return M,entwicklung

def u14_2_1():
    M = 25
    t = 20
    def P(x,t):
        D = 1
        d = 1
        return 1.0/np.sqrt(4*np.pi*D*t) * (np.exp(-(x-M)**2/(4.0*d*t)) - np.exp(- (x+M)**2/(4.0*d*t)))

    plt.plot([P(x,t) for x in range(0,2000)] ) 
    games = []
    laeufe = 2000
    for _ in range(laeufe):
        stand,_ = spielen(25, t)
        games.append(stand)
    plt.hist(games, normed=True)
    plt.show()
    
def u14_2_2():  
    M = 4
    D = 1 
    def m(t):
        return M/np.sqrt(2*D*t)
    
    def Q(t):
        m_ = m(t)
        return np.exp(-0.5*m_**2) * m_ * np.sqrt(2.0/np.pi) 
    
    t = 5000
    plt.plot(range(1,t), [Q(t) for t in range(1,t)])
    plt.show()
    
    # Spielen gegen die Bank
    for _ in range(300):
        _,entwicklung = spielen(M,np.inf)
        plt.plot(entwicklung)
    plt.title('Spielen gegen die Bank')
    plt.legend() 
    plt.show()
    
if __name__ == '__main__':
    u14_2_1()
    u14_2_2()
