#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Carlos Martín, Tran Tu

import numpy as np
import matplotlib.pyplot as plt
import random
from math import exp, sqrt, ceil
from multiprocessing import Pool
from itertools import product, groupby

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
    if random.uniform(0, 1) < min(1, chance):
        return True, dE

    return False, dE

def flip_spin(field, i, j):
    s = field[i,j]
    flip = +1.0 if s in [0.5, -1.5] else -1.0
    snew = s + flip

    field[i,j] = snew

# TODO: we should be able to calculate the diff here as well
def magnetisation(field):
    """Calculate the magnetisation of this field"""
    n = len(field)

    m = 0.0
    for i in range(n):
        for j in range(n):
            m += field[i,j]

    return m

# numpy.split() doesn't like it when we have lists as elements, so
# let's put them in a class
class FieldContainer:
    def __init__(self):
        self.fields = []

    def append(self, field):
        self.fields.append(field)
    

def simulate(n, kBT, H, steps, equilibrium, snapshot=False):
    """Simulate a field of size n**2"""

    field = random_field(n)
    E = total_energy(field, H)

    Es = [] # history of the energy in the system
    Ms = [] # history or magnetisations
    snaps = FieldContainer()

    for step in range(steps):
        # Generate a random pair, to see which spin we should try to
        # flip
        i, j = np.random.random_integers(0, n-1, 2)

        accepted, dE = accept_flip(E, H, kBT, field, i, j)
        #print accepted, dE

        if accepted:
            E += dE
            flip_spin(field, i, j)

        # Take a snapshot when we start counting and when we stop
        if snapshot and step == equilibrium + 1:
            snaps.append(field.copy())
        if snapshot and step == steps -1:
            snaps.append(field.copy())
            
            
        if step > equilibrium:
            Es.append(E)
            Ms.append(magnetisation(field))

    return np.mean(Es), np.mean(Ms), np.mean(np.abs(Ms)), snaps

# Run our simulation, we get the variables via a tuple though Pool.map()
def run_simulation((kBT, H)):
    n = 10
    nsteps = 25000
    nequi = nsteps / 4 # the first fourth is to let it settle

    return simulate(n, kBT, H, nsteps, nequi, True)

if __name__ == '__main__':
    repeats = 50 # repeats per experiment
    kBTs = [1e-3, 1e-1, 1, 2, 5, 10, 100]
    Hs = [0, 0.1, 1]
    # Cartesian product of the temperatures and the magenic fields
    experiments = list(product(kBTs, Hs))

    experiments = reduce(list.__add__, [[x]*repeats for x in experiments])      
    print 'experiments', experiments
    
    # Each experiment (and each repetition) is independent, so we can
    # create a pool of workers and make use of all the machine's
    # cores.
    pool = Pool()
    answers = pool.map(run_simulation, experiments)

    # The answers contain several copies of each experiment, so we
    # need to split it up into groups to average them
    deflattened = np.split(np.array(answers), len(kBTs)*len(Hs))

    # deflattened now has a bunch of arrays, each with tuples with the
    # energy and magnetism, and we need to extract them

    Es = []
    Ms = []
    AMs = []
    snaps = []
    for ex in deflattened:
        Es.append(np.mean(ex[:,0]))
        Ms.append(np.mean(ex[:,1]))
        AMs.append(np.mean(ex[:,2]))
        snaps.append(ex[:,3])

    print 'snaps', snaps
    # Plot the energy
    for i, h in enumerate(Hs):
        plt.semilogx(kBTs, Es[i::3], label='<E>, H = %g' % h)

    plt.title('Energie des Systems')
    plt.xlabel('$k_B T$')
    plt.ylabel('Energy')
    plt.legend()

    plt.figure()
    # Plot the magnetisms
    for i, h in enumerate(Hs):
        plt.semilogx(kBTs, Ms[i::3], label='<M> H = %g' % h)
        plt.semilogx(kBTs, AMs[i::3], label='<|M|> H = %g' % h)

    plt.title('Magnetismus des Systems')
    plt.xlabel('$k_B T$')
    plt.ylabel('<M>')
    plt.legend()

    # snaps is a list of FieldContainers (pairs)
    for containers, h in zip(snaps, Hs):
        side = ceil(sqrt(len(containers)))
        plt.figure()
        for i in range(len(containers)):
            plt.subplot(side, side, i)
            plt.xticks([])
            plt.yticks([])
            plt.imshow(containers[i].fields[0], extent=(0, 10, 10, 0), interpolation='nearest', cmap=plt.cm.jet)

        plt.suptitle('Snapshot bei H %g, am Equilibrium' % h)
        plt.savefig('snaps/h%s_equi.png' % (h))
        plt.figure()
        for i in range(len(containers)):
            plt.subplot(side, side, i)
            plt.xticks([])
            plt.yticks([])
            plt.imshow(containers[i].fields[1], extent=(0, 10, 10, 0), interpolation='nearest', cmap=plt.cm.jet)

        
        plt.suptitle('Snapshot bei H %g, am Ende' % h)
        plt.savefig('snaps/h%s_ende.png' % (h))


    #plt.show()
