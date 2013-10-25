import math
import numpy

def mround ( x,N ) :
    if ( x==0 ) :
        return x
    return round ( x, int ( N - math.ceil ( math.log10 ( abs ( x )))))

def leibnitz(k, n):
    sum = 0.0
    rev_sum = 0.0
    single = numpy.single(0.0)
    
    # 0 bis k
    for i in range(0,k):
        sum = sum + mround(pow(-1,i) / (2.0*i+1),n)
    r = range(0,k)
    
    list.reverse(r) 
    
    # umgekerte Reihenfolge
    for i in r:
        rev_sum = rev_sum + mround(pow(-1,i) / (2.0*i+1),n)
    # single precision 
    for i in range(0,k):
        single = single + numpy.single(pow(-1,i) / (2.0*i+1))
        
    return (sum, rev_sum, single)    

def stagniert(n):
    return math.ceil(math.sqrt(1/n))

#Aufgabe 1.1.1 und 1.1.2
print "%g %g %g" % leibnitz(500,4)
#>>> 0.784842 0.784842 0.784898
print "%g %g %g" % leibnitz(500,6)
#>>> 0.784898 0.784898 0.784898

#Aufgabe 1.1.3
# n = 6
print stagniert(0.0000004)
#>>> k = 1582.0    

# n = 4
print stagniert(0.00004)
#>>> k = 159.0
