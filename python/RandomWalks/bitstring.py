#!/usr/bin/env python

import numpy as np
import random
import itertools

############################################################
# GA (bitstring) stuff
############################################################

def hamming_distance(x, y):
    return np.sum(x != y)

def bitstring_per_ind_mutation(x):
    i = random.randrange(len(x))
    x[i] = not x[i]
    return x
    
def make_bitstring_per_gene_mutation(pmut):
    def m(x):
        for i in range(len(x)):
            if random.random() < pmut:
                x[i] = not x[i]
        return x
    return m

def per_ind_sd(n):
    """Calculate SD for a per-individual mutation given the bitstring
    length n

    FIXME one of these is not right
    """
    n = float(n)
    N = 2**n
    mu = 1.0/N
    nneighbours = n
    
    # we have n neighbours, each has p = 1/n
    A = ((1.0/(2**n) - 1.0/n)**2) * n
    # we have 2**n - n non-neighbours, each has p = 0
    B = ((1.0/(2**n) - 0)**2) * (2**n - n)
    
    var = (A + B) / (2**n)
    
    var2 = n / (2 ** (3*n)) - 1.0 / (2**(2*n-1)) + 1.0 / (n * (2**n)) - n / (2**(2*n))

    
    return math.sqrt(var), math.sqrt(var2)

def generate_bitstring_tm_row0(length, pmut=None):
    row0 = np.zeros(2**length)
    i = 0
    indi = [0] * length
    for j, indj in enumerate(itertools.product(*[(0, 1) for y in range(length)])):
        indj = np.array(indj, dtype='bool')
        h = hamming_distance(indi, indj)
        if pmut is None:
            if h == 1:
                row0[j] = 1.0 / length # there are length inds at hamming distance 1
                # else leave it at zero
        else:
            row0[j] = (pmut ** h) * ((1.0 - pmut) ** (length - h))
    return row0
    
def generate_ga_tm(length, pmut=None):
    """For a bitstring (genetic algorithm) representation of a given
    length, generate a transition matrix with the mutation probability
    pmut. Also generate the Hamming distances. If pmut=None (default),
    exactly one bitflip is performed per individual, rather than using
    a per-gene mutation probability."""

    tm = np.zeros((2**length, 2**length))
    hm = np.zeros((2**length, 2**length))
    for i, indi in enumerate(itertools.product(*[(0, 1) for x in range(length)])):
        indi = np.array(indi, dtype='bool')
        for j, indj in enumerate(itertools.product(*[(0, 1) for y in range(length)])):
            indj = np.array(indj, dtype='bool')
            h = hamming_distance(indi, indj)
            hm[i][j] = h
            if pmut is None:
                if h == 1:
                    tm[i][j] = 1.0 / length # there are length inds at hamming distance 1
                    # else leave it at zero
            else:
                tm[i][j] = (pmut ** h) * ((1.0 - pmut) ** (length - h))
    return tm, hm

def nCk(n, k):
    """n-choose-k"""
    return scipy.misc.comb(n, k, True)

def Krovi_Brun_bitstring_MFPT(n, d):
    """Calculate the MFPT between two bitstrings. n is the bitstring
    length, and d is the Hamming distance between two individuals.
    Formula is given in
    http://math.stackexchange.com/questions/28179/logic-question-ant-walking-a-cube/28188#28188,
    also a variation is given in 'Hitting time for quantum walks on
    the hypercube', Hari Krovi and Todd A. Brun, PHYSICAL REVIEW A 73,
    032341 2006. This is not needed for anything else, just as a check
    that the formula gives same answers as the general MFPT formula
    (see ergodic.py) -- it does."""
    return sum(
        sum(nCk(n, k) for k in range(m+1)) / float(nCk(n-1, m))
        for m in range(n-d, n))

def onemax_fitvals(length):
    return [np.sum(ind) for ind in
            itertools.product(*[(0, 1) for x in range(length)])]

def ga_tm_wrapper(dirname, pmut=None):
    # dirname should be <dir>/ga_length_6, for example
    length = int(dirname.strip("/").split("_")[2])
    tm, hm = generate_ga_tm(length, pmut)
    outfilename = dirname + "/TP.dat"
    np.savetxt(outfilename, tm)
    outfilename = dirname + "/Hamming.dat"
    np.savetxt(outfilename, hm)

