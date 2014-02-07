#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín, Tran Tu

import numpy as np
import random

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
    #vals = (x for x in np.random.random_integers(0, 1, max_conns))

    # Generate the n-by-n edge matrix
    M = np.empty((n,n), dtype=bool)
    M.fill(False)
    # and fill it with the random values from before
    for i in range(1, n):
        for j in range(i):
            M[j,i] = next(vals)

    # We now need to mirror across the diagonal
    for i in range(1, n):
        for j in range(i):
            M[i,j] = M[j,i]

    # transform what we return into a boolean array
    return M

if __name__ == '__main__':
    M = random_network(5, 0.5)
    print 'M\n', M
