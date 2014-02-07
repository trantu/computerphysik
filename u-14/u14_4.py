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

    # Generate the n-by-n edge matrix
    M = np.zeros((n,n), dtype=int)
    # and fill it randomly with 1 or 0
    for j in range(1, n):
        for i in range(j):
            if random.random() <= p:
                M[i,j] = 1

    # We now need to mirror across the diagonal. This makes it easier
    # to deal with edges, as we don't need to worry which is the entry
    # that contains the true value.
    M = M + M.transpose()

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
        return c

if __name__ == '__main__':
    N = 1000
    p = 0.02

    random.seed()
    
    dist = [0]*N
    # Count connections 10 times
    pool = Pool()
    lists = pool.map(one_random, [(N, p)]*10)
    counts = [0]*N
    # and average how many neighbours each node had
    for i in range(N):
        counts[i] = np.mean([lists[j][i] for j in range(10)])
        
    norm, bins, _ = plt.hist(counts, normed=True, label='$P(k_i)$')

    pois = poisson(N, p)
    plt.hist(pois, normed=True, label='Poisson', alpha=0.5)

    plt.legend()
    plt.show()
