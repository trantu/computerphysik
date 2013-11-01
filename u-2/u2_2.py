'''
Created on Nov 1, 2013

@author: tran
'''
import math
import numpy

def diff_quot(xs, ys):
    res = [((ys[i+1]-ys[i-1])/(xs[i+1]-xs[i-1]),-(math.sin(xs[i]))) for i in range(1,len(ys)-1)]
    return res

def k():
    delta_x = [2*math.pi/2**k for k in range(100)]
    #print delta_x
    counter = 0
    for i in delta_x:
        #print i
        xs = numpy.arange(0.0,100,i)
        ys = map(math.cos, xs)
        abs_fails = [abs(a-b) for (a,b) in diff_quot(xs, ys)]
        print abs_fails
        for j in abs_fails:
            if j >= 0.01:
                print "Fehler = %g , Delta_x = %g" %(j,i)
                return counter
        counter += 1           
    return counter
#print res[counter]
print "k = %i" % k()

'''
Fehler = 0.36338 , Delta_x = 1.5708
k = 2
'''

