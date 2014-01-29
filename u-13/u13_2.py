#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Carlos Martín, Tran Tu

import numpy as np
import random
from math import exp

random.seed(5)

# Possible spins
spins = [+0.5, -0.5, +1.5, -1.5]

# Sie of each side
N = 10

def random_field(n):
    field = np.zeros((N, N), dtype=float)
    # Fill it with random spins from our list of possible ones
    for i in range(n):
        field[i] = [spins[int(random.uniform(0, 3))] for _ in range(N)]

    return field

def total_energy(field, H):
    """Calculate the total energy in the given field.

    field the magnetic field
    H strenght of the field
    We assume that the field is square"""

    n = len(field)

    lhs = rhs = 0.0
    # Since we're working with J_ij=0 for non-neighbours, we can
    # simply ignore them in the left-hand sum. That is, we only
    # calculate s_i*s_j for the four immediate neighbours.
    #
    # The formula uses only one axis, so we use that index as the
    # linearised version of our field.
    for i in range(n):
        for j in range(n):
            higher = field[(i-1)%n,j]
            lower  = field[(i+1)%n,j]
            left   = field[i,(j-1)%n]
            right  = field[i,(j+1)%n]
            # This is the same as multiplying each pair and adding them up
            lhs += field[i,j] * (higher + lower + left + right)
            rhs += field[i,j]

    return -0.5 * lhs - H * rhs

def energy_change(field, i, j, h):
    """Calculate the change in energy if we flip the spin at position (i,j)"""

    n = len(field)
    s = field[i,j]
    flip = +1.0 if s in [0.5, -1.5] else -1.0
    snew = s + flip

    higher = field[(i-1)%n,j]
    lower  = field[(i+1)%n,j]
    left   = field[i,(j-1)%n]
    right  = field[i,(j+1)%n]

    e_old = -s * (higher + lower + left + right) - h * s
    e_new = -snew * (higher + lower + left + right) - h * snew

    return e_new - e_old

def accept_flip(E, H, kBT, field, i, j,):
    dE = energy_change(field, i, j, H)

    # we always accept a flip that loses energy
    if dE < 0:
        return True, dE

    # if it doesn't lose energy, there's still a chance to accept it
    E2 = E + dE
    chance = exp((E - E2) / kBT)
    if random.uniform(0, 1) < chance:
        return True, dE

    return False, dE

def flip_spin(field, i, j):
    s = field[i,j]
    flip = +1.0 if s in [0.5, -1.5] else -1.0
    snew = s + flip

    field[i,j] = snew

def simulate(n, H, kBT, steps=100):
    """Simulate a field of size n**2"""

    field = random_field(n)
    E = total_energy(field, H)

    for step in range(steps):
        # Generate a random pair, to see which spin we should try to
        # flip
        i, j = np.random.random_integers(0, n-1, 2)

        accepted, dE = accept_flip(E, H, kBT, field, i, j)
        print accepted, dE

        if accepted:
            E += dE
            flip_spin(field, i, j)

if __name__ == '__main__':
    simulate(10, 1, 10e-3)
