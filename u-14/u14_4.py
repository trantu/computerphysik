#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín, Tran Tu

import numpy as np
import matplotlib.pyplot as plt
import random
from math import exp
from math import factorial as fac
from collections import Counter
from multiprocessing import Pool
import time

def random_network(n, p):
    """Create a randomly connected n-node network.

    n: size of the network
    p: probability that two nodes are connected

    Returns a boolean array of edges"""

    rng = random.Random()

    # first, we need n(n-1)/2 random values {0,1} to represent the
    # connections We make 'vals' a generator so we can keep reading
    # data from it, as we don't care about position.
    max_conns = n*(n-1) / 2
    vals = (x < p for x in [rng.uniform(0, 1) for _ in range(max_conns)])

    # Generate the n-by-n edge matrix
    M = np.empty((n,n), dtype=bool)
    M.fill(False)
    # and fill it with the random values from before
    for i in range(1, n):
        for j in range(i):
            M[j,i] = next(vals)

    # We now need to mirror across the diagonal. This makes it easier
    # to deal with edges, as we don't need to worry which is the entry
    # that contains the true value.
    for i in range(1, n):
        for j in range(i):
            M[i,j] = M[j,i]

    # transform what we return into a boolean array
    return M

def count_edges(M):
    """Return the amount of edges in the network M"""

    # As the values are repeated, get only a triangle
    return np.count_nonzero(np.tril(M))

def count_conns(M, i):
    """Count the amount of connections which node 'i' has"""

    # we can simply count non-false for the row/column
    return np.count_nonzero(M[i])

def poisson(N, p):
    lam = p * (N - 1) # expected value <k>

    # NumPy helpfully provides us with a Poisson function which
    # doesn't complain about longs being too big.
    return np.random.poisson(lam=lam, size=N)

def one_random((n, p)):
        M = random_network(N, p)
        c = [count_conns(M, i) for i in range(N)]
        return Counter(c).items()

if __name__ == '__main__':
    N = 1000
    p = 0.02

    random.seed()
    
    dist = [[0]*N]*10
    # Count connections 10 times
    pool = Pool()
    for n, items in enumerate(pool.map(one_random, [(N, p)]*10)):
        print(items)
        print(dist[9])
        for i, reps in items:
            dist[n][i] += reps

    # Now that we have the ten repetitions, we need to find the
    # average amount of nodes with a particular amount of neighbours.

    dist = np.matrix(dist, dtype=int)
    avg_conns = []
    for i in range(N):
        avg_conns.append(np.mean(dist[:,i]))

    plt.hist(avg_conns, bins=100, normed=True, label="Real")

    pois = poisson(N, p)
    plt.hist(pois, bins=100, normed=True, label='Poisson')

    plt.legend()
    plt.show()
