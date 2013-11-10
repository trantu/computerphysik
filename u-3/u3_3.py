#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto, Tu Tran

from numpy import matrix
from numpy.linalg import solve, inv
from math import sqrt

def frobenius(A):
    """Calculate the Frobenius norm of `A`. We assume that the matrix is
    square and is of numpy.matrix type.

    """

    n, m = A.shape
    # square each element, returning them in a flattened list
    squared = [abs(A[i,j])**2 for i in range(n) for j in range(m)]
    return sqrt(sum(squared))

def show_solution(A, b, bp, ls):
# "real" solution
    x = solve(A, b)

    for l in ls:
        bc = l * bp # current b matrix
        db = abs(bc - b)

        print "With %f: " % l
        print "\t║𝚫b║₂ = %f" % frobenius(db)

        # solve it
        xc = solve(A, l*bp)
        dx = abs(xc - x)

        # Calcualte the errors
        Ainv = inv(A)
        cond = frobenius(A) * frobenius(Ainv)
        abserr = frobenius(Ainv) * frobenius(db)
        relerr = cond * (frobenius(db) / frobenius(b))

        print "\t Expected max error %f" % abserr
        print "\t║𝚫x║₂ = %f" % frobenius(dx)
        print "\t Expected rel error %f" % relerr
        print "\t║𝚫x║₂/║x║₂ = %f" % (frobenius(dx) / frobenius(x))
        print ""

if __name__ == '__main__':
    ## Set up the data per the assignment text
    A1 = matrix([[1,       -2, -1],
                 [36000,    2,  0],
                 [-2,    1400,  1]])

    b1 = matrix([[3.0], [72002.0], [1399.0]])
    b1p = matrix([[1.1  * b1[0,0]],
                  [0.9  * b1[1,0]],
                  [1.05 * b1[2,0]]])

    A2 = matrix([[1.0,   1.0/2, 1.0/3],
                 [1.0/2, 1.0/3, 1.0/4],
                 [1.0/3, 1.0/4, 1.0/5]])

    rowsum = A1.sum(1, dtype='float').tolist()

    b2 = matrix([rowsum[0],
                 rowsum[1],
                 rowsum[2]])

    b2p = matrix([[1.1  * b2[0,0]],
                  [0.9  * b2[1,0]],
                  [1.05 * b2[2,0]]])

    
    ls = [1.0, 0.1, 5.0]

    print "A1"
    show_solution(A1, b1, b1p, ls)

    print "A2"
    show_solution(A2, b2, b2p, ls)
