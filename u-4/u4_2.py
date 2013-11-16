#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from mpl_toolkits . mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

from numpy import vectorize, matrix, allclose
from numpy.linalg import solve
from math import sqrt
from itertools import starmap
import sys

def phi(x, y):
    return x**4 - x**2 + y**4 - 0.2 * y**3 - y**2 + 0.2 * x * y**3

# $-\nable \phi$
def F(x, y):
    return matrix([[4 * x**3 - 2 * x + 0.2 * y**3],
                   [4 * y**3 - 0.6 * y**2 + 0.6 * x * y**2 - 2 * y]], dtype=float)

# Jacobi matrix for phi2
def DF(x, y):
    return matrix([[12 * x**2 - 2, 0.6 * y**2],
                   [0.6 * y**2, 12 * y**2 - 1.2*y + 1.2*x*y - 2]])

## Plot $\phi(x,y)$
x = np.arange ( -1.0 ,1.01 ,0.05)
y = np.arange ( -1.0 ,1.01 ,0.05)
x, y = np.meshgrid (x , y )

#a = Axes3D(plt.figure())
#a.plot_wireframe(x, y, np.vectorize(phi)(x, y))
#plt.show()

def norm2(A):
    n, m = A.shape
    nums = [A[i,j]**2 for i in range(n) for j in range(m)]
    return sqrt(sum(nums))

## Newton method for linear equations
def newton_next(f, D, vec):
    x = vec[0,0]
    y = vec[1,0]

    return solve(D(x,y), -f(x, y)) + vec

def newton(x, y, d):
    point = matrix([[x],
                    [y]], dtype=float)

    while not ((norm2(point) < d) or (norm2(F(point[0,0], point[1,0])) < d)):
        last_point = point
        point = newton_next(F, DF, last_point)
        #print point[0,0], point[1,0], norm2(point), norm2(F(point[0,0], point[1,0])), F(point[0,0], point[1,0])

        
    return point


# Calculate the guesses
guesses = [(1.0,  1.0),
           (1.0,  -1.0),
           (-1.0, 1.0),
           (-1.0, -1.0),
           (0.01,   0.01)]

delta = 0.001
extrema = [newton(x, y, delta) for (x, y) in guesses]

min1 = extrema[0]
min2 = extrema[1]
min3 = extrema[2]
min4 = extrema[3]
maxi = extrema[4]

def Zuordnung(x, y, f=newton):
    p = f(x, y, delta)
    if norm2(p - min1) < delta:
        return 1
    if norm2(p - min2) < delta:
        return 2
    if norm2(p - min3) < delta:
        return 3
    if norm2(p - min4) < delta:
        return 4
    if norm2(p - maxi) < delta:
        return 5

    return 0

## Get all points in the grid and plot them against where they end up
x = np.arange ( -1.0 ,1.01 ,0.02)
y = np.arange ( -1.0 ,1.01 ,0.02)
x, y = np.meshgrid (x , y )
Zuordnungvec = vectorize(Zuordnung)
#plt.imshow(Zuordnungvec(x, y), extent = [-1, 1, 1, -1])
#plt.colorbar()
#plt.show()

def iteration_step(xy, e, F):
    x = xy[0,0]
    y = xy[1,0]

    return xy + e * F(x, y)

def iteration(x, y, f, e, d):
    point = matrix([[x], [y]], dtype=float)

    while not norm2(f(point[0,0], point[1,0])) < d:
        point = iteration_step(point, e, f)

    return point

## \epsilon = -D[f(x0)]^-1 ist ähnlich zum vereinfachtes Newton.
Dinv = -np.linalg.inv(DF(-1.0, 1.0))
def find_min(x, y, d):
    sys.stdout.write("%f, %f           \r" % (x, y))
    return iteration(x, y, F, Dinv, d)

def Zuordnung2(x, y):
    return Zuordnung(x, y, f=find_min)

print find_min(1.0, 1.0, delta)

x = np.arange ( -1.0 ,1.01 ,0.02)
y = np.arange ( -1.0 ,1.01 ,0.02)
x, y = np.meshgrid (x , y )
Zuordnungvec = vectorize(Zuordnung2)
plt.imshow(Zuordnungvec(x, y), extent = [-1, 1, 1, -1])
plt.colorbar()
plt.show()
