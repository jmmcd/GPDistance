#!/usr/bin/env python

"""This module provides functionality implementing an absorbing Markov
chain. It can be used to calculate the expected length of a random
walk between any two points in a space, given the transition matrix on
the space. 

"""

import numpy as np
import random
import sys

# this import is from
# http://www.pysal.org/library/spatial_dynamics/ergodic.html: it
# implements the first-mean-passage-time algorithm, as described in
# "Introduction to Probability Theory with Computing" by J. Laurie
# Snell (Prentice-Hall, 1975).
import ergodic 

## normalise an array row-by-row, ie make each row sum to 1. This
## might be needed for matrices created using a hill-climbing
## adjustment to transition probabilies: in such cases, the transition
## probabilities are adjusted to take account of fitness, so "bad"
## transitions (to worse fitness) are made less likely to be accepted,
## meaning that each row no longer sums to 1. By renormalising we get
## the true probability after (if necessary) multiple rounds of
## rejection and finally one acceptance.
def normalise(d):
    for i in range(len(d)):
        sumval = np.sum(d[i])
        d[i] *= 1.0 / sumval
    return d

def make_random_matrix(n):
    """Make a random transition matrix on n points. For testing
    only.

    """
    tm = np.zeros((n, n))
    for i in range(n):
        # generate a random vector of out-probabilities
        vec = np.zeros((n, 1))
        for j in range(n):
            vec[j] = random.random()
        # normalize
        vec *= 1.0 / sum(vec)
        for j in range(n):
            tm[i, j] = vec[j]
    return tm

def make_absorbing(tm, dest):
    """Given a transition matrix, disallow transitions away from the
    given destination -- ie make it an absorbing matrix.

    """
    for j in range(len(tm)):
        tm[dest, j] = 0.0
    tm[dest, dest] = 1.0
    return tm

def read_transition_matrix(filename):
    """Read a transition matrix from a file and return. The matrix
    will have been written in the right format by some Java code.

    """
    d = np.genfromtxt(filename)
    # check that each row sums to 1, since each row is the
    # out-probabilities from a single individual. Allow a small margin
    # of error.
    epsilon = 0.000001
    for i, v in enumerate(d):
        assert(1.0 - epsilon < sum(v) < 1.0 + epsilon)
    return d
    
def run_simulation(tm, src, dest, conf, max_iters):
    """Given a transition matrix, a source and destination index, and
    a confidence, calculate the number of iterations required for the
    src to transition to the destination with the given confidence. If
    max_iters is exceeded, return max_iters.

    """
    n = len(tm)
    initial = np.zeros((1, n))
    initial[0, src] = 1.0

    iters = 0
    finished = False
    s = dot(initial, tm)
    while not finished:
        # print("iters = {0}, dest = {1}".format(iters, s[0, dest]))
        if iters >= max_iters:
            finished = True
        if s[0, dest] > conf:
            finished = True
        s = dot(s, tm)
        iters += 1
    return iters

def run_many_simulations(filename, conf):
    """Given a transition matrix, calculate the number of iterations
    required to move between all sources and destinations, and save as
    a new matrix.

    """

    tm = read_transition_matrix(filename)
    n = len(tm)
    max_iters = 5 * n
    retval = np.zeros((n, n))
    for dest in range(n):
        tm = read_transition_matrix(filename)
        tm = make_absorbing(tm, dest)
        for src in range(n):
            print("n = {0}, src = {1}, dest = {2}".format(n, src, dest))
            iters = run_simulation(tm, src, dest, conf, max_iters)
            retval[src, dest] = iters
    np.savetxt("depth_2_rw_conf_0.5.dat", retval)
    

def run_test():
    # standard parameters
    conf = 0.5

    use_random = False
    if use_random:
        # for testing only
        n = 1000
        tm = make_random_matrix(n)
    else:
        tm = read_transition_matrix("depth_2_tm.dat")
        n = len(tm)

    max_iters = 5 * n

    # arbitrary
    src = random.randrange(n)
    dest = random.randrange(n)

    tm = make_absorbing(tm, dest)

    print("n = {0}, src = {1}, dest = {2}".format(n, src, dest))
    print(tm)
    run(tm, src, dest, conf, max_iters)

def get_fmpt(x):
    x = np.matrix(x) # NB! has to be a true matrix or it breaks somehow...
    return ergodic.fmpt(x)
    
def read_and_get_fmpt(infilename, outfilename):
    x = read_transition_matrix(infilename)
    f = get_fmpt(x)
    np.savetxt(outfilename, f)

def generate_random_tm(n):
    d = np.zeros((n, n), dtype=float)

    for i in range(len(d)):
        v = np.random.random(len(d))
        mv = np.sum(v)
        v *= 1.0 / mv
        d[i] = v
    d = np.matrix(d)
    return d

def test_matrix_size(n):
    """Test how big the tm can be before get_fmpt becomes
    slow. n = 4000 is fine, n = 10000 starts paging out (at least 30
    minutes).

    """
    d = generate_random_tm(n)
    fmpt = get_fmpt(d)
    print("min", np.min(fmpt))
    print("max", np.max(fmpt))

def invert_probabilities(adj):
    assert(np.min(adj) > 0.0)
    return -np.log(adj)

def deinvert_probabilities(adj):
    return np.exp(-adj)

def test_floyd_warshall_random_data(n):
    adj = generate_random_tm(n)
    return floyd_warshall(adj)
    
def floyd_warshall(adj):
    """Finds the shortest path between all pairs of nodes. For this to
    be useful, need to invert the transition matrix probabilities p
    somehow, eg using -log(p), so that low probabilities cause high
    edge costs. Can then invert the shortest path values s again using
    e^(-s) to recover a probability. n = 1000 works fine, n = 4000
    starts paging out (at least 10 minutes).

    """
    adj = invert_probabilities(adj)
    # from
    # http://www.depthfirstsearch.net/blog/2009/12/03/computing-all-shortest-paths-in-python/
    n = len(adj)
    for k in range(n):
        adj = np.minimum(adj,
                         np.tile(adj[:, k].reshape(-1, 1), (1, n)) +
                         np.tile(adj[k, :], (n, 1)))
    adj = deinvert_probabilities(adj)
    return adj


def read_and_get_fmpt_hpp(codename):
    d = read_transition_matrix(codename + "1STP.dat")
    f = get_fmpt(d)
    outfilename = codename + "FMPT.dat"
    np.savetxt(outfilename, f)

    h = floyd_warshall(d)
    outfilename = codename + "HPP.dat"
    np.savetxt(outfilename, h)
    
if __name__ == "__main__":
    codename = sys.argv[1]
    read_and_get_fmpt_hpp(codename)