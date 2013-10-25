#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Carlos Martín Nieto
# Tu Tran

def find_e(a):
    e = 1.0

    while a != a + e:
        e /= 2

    return e

# Store the $\alpha$ and its $\epsilon$ in the tuple directly so we
# can simply print it
es = [(a, find_e(a)) for a in [1.0, 1.1e-5, 1e-30, 5e12]]

for val in es:
    print "%g\t%g" % val
