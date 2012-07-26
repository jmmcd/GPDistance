#!/usr/bin/env python

"""This module provides functionality implementing an absorbing Markov
chain. It can be used to calculate the expected length of a random
walk between any two points in a space, given the transition matrix on
the space. Also the lowest-cost path (treating low transition
probabilities as high edge traversal costs) and the shortest path (in
number of steps, disregarding probabilities).

"""

import numpy as np
import random
import sys
import itertools

# this import is from
# http://www.pysal.org/library/spatial_dynamics/ergodic.html: it
# implements the first-mean-passage-time algorithm, as described in
# Kemeny, John, G. and J. Laurie Snell (1976) Finite Markov
# Chains. Springer-Verlag, Berlin.
import ergodic

def normalise_by_row(d):
    """Normalise an array row-by-row, ie make each row sum to 1. This
    is useful when randomly-generating transition matrices, and for
    matrices created using a hill-climbing adjustment to transition
    probabilies: in such cases, the transition probabilities are
    adjusted to take account of fitness, so "bad" transitions (to
    worse fitness) are made less likely to be accepted, meaning that
    each row no longer sums to 1. By renormalising we get the true
    probability after (if necessary) multiple rounds of rejection and
    finally one acceptance.

    """
    for i in range(len(d)):
        sumval = np.sum(d[i])
        d[i] *= 1.0 / sumval
    return d

def make_random_matrix(n):
    """Make a random transition matrix on n points."""
    tm = np.random.random((n, n))
    return normalise_by_row(tm)

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
    """Run a simulation with test parameters."""
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
    run_simulation(tm, src, dest, conf, max_iters)

def set_self_transition_zero(x):
    """Set cost/length of self-transition to zero."""
    # Was using this method:
    #
    # x *= (np.ones_like(x) - np.eye(len(x)))
    #
    # ... but it fails if infinities are involved. The following
    # method works ok.
    for i in range(len(x)):
        x[i][i] = 0.0
        

def get_fmpt(x):
    """Calculate first-mean-passage time of a given transition
    matrix. Set self-transitions to zero.

    """
    # NB! The ergodic code expects a matrix, not a numpy array. Breaks
    # otherwise.
    x = np.matrix(x)
    x = np.array(ergodic.fmpt(x))
    set_self_transition_zero(x)
    return x
    
def test_matrix_size(n):
    """Test how big the tm can be before get_fmpt becomes
    slow. n = 4000 is fine, n = 10000 starts paging out (at least 30
    minutes).

    """
    d = make_random_matrix(n)
    fmpt = get_fmpt(d)
    print("min", np.min(fmpt))
    print("max", np.max(fmpt))

def invert_probabilities(adj):
    """Convert a probability into an edge traversal cost. It's ok to
    ignore runtime floating point problems here, because they're only
    due to zero-probabilities, which result in infinite edge-traversal
    costs, which is what we want. Restore the "raise" behaviour
    afterward.

    """
    np.seterr(all='warn')
    retval = -np.log(adj)
    np.seterr(all='raise')
    return retval

def deinvert_probabilities(adj):
    """Convert an edge traversal cost into a probability."""
    return np.exp(-adj)

def test_floyd_warshall_random_data(n):
    """Test."""
    adj = make_random_matrix(n)
    return floyd_warshall_probabilities(adj)

def floyd_warshall_probabilities(adj):
    """For this to be useful, need to invert the transition matrix
    probabilities p somehow, so that low probabilities cause high edge
    traversal costs. See invert_ and deinvert_probabilities.

    """
    x = invert_probabilities(adj)
    x = floyd_warshall(x)
    set_self_transition_zero(x)
    return x
    
def floyd_warshall(adj):
    """Finds the shortest path between all pairs of nodes. For this to
    be useful, the edge weights have to have the right semantics: a
    small transition probability implies a large edge weight
    (cost). So if your edge weights are probabilities, use
    floyd_warshall_probabilities() instead: it performs the
    conversion.

    Concerning matrix size: n = 1000 works fine, n = 4000 starts
    paging out (at least 10 minutes) on my machine.

    """
    # from
    # http://www.depthfirstsearch.net/blog/2009/12/03/computing-all-shortest-paths-in-python/
    n = len(adj)
    for k in range(n):
        adj = np.minimum(adj, np.add.outer(adj[:,k],adj[k,:]))
    return adj


def floyd_warshall_nsteps(adj):
    """Disregard the transition probabilities, other than to see
    whether an edge traversal is allowed or not. Calculate the number
    of steps required to get from each point to each other.

    """
    x = discretize_probabilities(adj)
    x = floyd_warshall(x)
    set_self_transition_zero(x)
    return x
    
def discretize_probabilities(d):
    """Set the edge cost to 1 if there is a nonzero probability, and
    to infinity if there is a zero probability. 

    """
    retval = np.ones_like(d, dtype=float)
    inf = np.infty
    for i in range(len(d)):
        for j in range(len(d)):
            if not d[i, j] > 0.0:
                retval[i, j] = inf
    return retval

def get_dtp(t):
    x = invert_probabilities(t)
    set_self_transition_zero(x)
    return x
    
def read_and_get_dtp_fmpt_sp_steps(codename):
    t = read_transition_matrix(codename + "/TP.dat")

    # This gets D_TP, which is just the transition probability inverted
    d = get_dtp(t)
    outfilename = codename + "/D_TP.dat"
    np.savetxt(outfilename, d)
    
    # This gets the first mean passage time, ie the expected length of
    # a random walk.
    f = get_fmpt(t)
    outfilename = codename + "/FMPT.dat"
    np.savetxt(outfilename, f)

    # This gets the cost of the shortest path between pairs. The cost
    # of an edge is the negative log of its probability.
    h = floyd_warshall_probabilities(t)
    outfilename = codename + "/SP.dat"
    np.savetxt(outfilename, h)

    # this gets the minimum number of steps required to go between
    # pairs, disregarding probabilities. Only interesting if some
    # edges are absent (ie edge probability is zero).
    p = floyd_warshall_nsteps(t)
    outfilename = codename + "/STEPS.dat"
    np.savetxt(outfilename, p)

def hamming_distance(x, y):
    return np.sum(x != y)
    
def generate_ga_tm(codename, pmut=None):
    """For a bitstring (genetic algorithm) representation of a given
    length, generate a transition matrix with the mutation probability
    pmut. Also generate the Hamming distances. If pmut=None (default),
    exactly one bitflip is performed per individual, rather than using
    a per-gene mutation probability.

    """
    length = int(codename.split("_")[2])
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
    outfilename = codename + "/TP.dat"
    np.savetxt(outfilename, tm)
    outfilename = codename + "/Hamming.dat"
    np.savetxt(outfilename, hm)

    
if __name__ == "__main__":
    codename = sys.argv[1]
    if "per_ind" in codename:
        generate_ga_tm(codename)
    else:
        generate_ga_tm(codename, 0.1)
    read_and_get_dtp_fmpt_sp_steps(codename)
