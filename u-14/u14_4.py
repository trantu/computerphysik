#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín, Tran Tu

import numpy as np
import matplotlib.pyplot as plt
import random
from math import exp
from math import factorial as fac
from collections import deque
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

    # we can simply count non-zero for the row/column
    return np.count_nonzero(M[i])

def poisson(N, p):
    lam = p * (N - 1) # expected value <k>

    # NumPy helpfully provides us with a Poisson function which
    # doesn't complain about longs being too big.
    return np.random.poisson(lam=lam, size=N)

def one_random((n, p)):
        M = random_network(n, p)
        c = [count_conns(M, i) for i in range(n)]
        return c

def neighbours(M, i):
    """Return a list of *i*'s neighbours"""

    ret = []
    for j, val in enumerate(M[i,:]):
        if val == 1:
            ret.append(j)

    return ret

def distances(M, i):
    # Keep a list of the nodes we've visited so we can avoid visiting
    # them again
    visited = [False] * len(M)
    visited[i] = True # we've visited ourselves
    
    # Distance from us to a particular node
    distance = [0] * len(M)
    # Nodes we should visit, start with ourselves
    q = deque([i])

    while len(q) > 0:
        k = q.popleft()
        neigh = neighbours(M, k)

        # Go through each of the node's neighbours which we've yet to
        # visit
        for n in neigh:
            if visited[n]:
                continue

            # Mark as visited and set the distance to that node as the
            # distance to the current node plus one
            visited[n] = True
            distance[n] = distance[k] + 1
            q.append(n)

    return distance

def u14_4_2():
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

def one_avg_distance((n, p)):
    M = random_network(n, p)
    return [distances(M, j) for j in range(n)]

def average_distance(n, pN):
    """Return the average distance, over five runs, of an n-node network"""
    p = float(pN) / n

    pool = Pool()
    dist = pool.map(one_avg_distance, [(n,p)]*5)
    return np.mean(dist)

if __name__ == '__main__':
    # u14_4_2()

    Ns = range(60, 161, 20)


    for pN in [10, 20]:
        res = []
        for N in Ns:
            res.append(average_distance(N, pN))
        plt.plot(Ns, res, label=r'$\langle l \rangle_N$, $pN = %d$' % (pN))

    plt.semilogy(Ns, Ns, label=r'$\log N$')
    plt.legend(loc='upper left')
    plt.show()
