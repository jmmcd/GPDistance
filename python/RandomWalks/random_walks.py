#!/usr/bin/env python

"""This module provides functionality implementing an absorbing Markov
chain. It can be used to calculate the expected length of a random
walk between any two points in a space, given the transition matrix on
the space. Also the lowest-cost path (treating low transition
probabilities as high edge traversal costs) and the shortest path (in
number of steps, disregarding probabilities)."""

import numpy as np
import scipy.stats
import random
import sys
import os
import itertools

# this import is from
# http://www.pysal.org/library/spatial_dynamics/ergodic.html: it
# implements the mean-first-passage-time algorithm, also known as the
# first-mean-passage-time, as described in Kemeny, John, G. and
# J. Laurie Snell (1976) Finite Markov Chains. Springer-Verlag,
# Berlin.
import ergodic

def analyse_random_walk(dirname):
    """Java code will write out a list of sampled lengths of random
    walks between nodes i and j. A single file, each line containing
    the tree i, tree j, then the list of samples. Analyse the
    basics."""
    
    n = 20
    f = open(dirname + "/MFPT_random_walking_samples.dat")
    x_mean = np.zeros((n, n))
    x_var = np.zeros((n, n))
    x_len = np.zeros((n, n), dtype=int)
    i = 0
    j = 0
    for line in f:
        t0, t1, samples = line.split(":")
        samples = np.array(map(int, samples.strip().split(" ")))
        x_mean[i, j] = np.mean(samples)
        x_var[i, j] = np.std(samples)
        x_len[i, j] = len(samples)
        i += 1
        if i == n:
            i = 0
            j += 1
    print(x_mean)
    print(x_var)
    print(x_len)

def estimate_MFPT_with_supernode(dirname):
    """Given a directory, go in there and get all the files under
    TP_supernode_estimates. For each, run the algorithm to get mfpt,
    and from the result extract the (0, 1) and (1, 0) values. Find out
    which trees they correspond to. Then correlate between all these
    values and the values calculated by exact MFPT given the complete
    TP matrix."""
    n = 50
    all_trees = open(dirname + "/all_trees.dat").read().strip().split("\n")
    estimate = np.zeros(2 * n)
    exact = np.genfromtxt(dirname + "/MFPT.dat")
    exact_extract = np.zeros(2 * n)
    ted = np.genfromtxt(dirname + "/TED.dat")
    ted_extract = np.zeros(2 * n)
    for i in range(n):
        d = read_transition_matrix(dirname + "/TP_supernode_estimates/"
                                   + str(i) + "_TP_estimates.dat")
        m = get_mfpt(d)
        t, s = open(dirname + "/TP_supernode_estimates/"
                    + str(i) + "_trees.dat").read().strip().split("\n")
        ti = all_trees.index(t)
        si = all_trees.index(s)
        estimate[2*i] = m[0, 1]
        estimate[2*i+1] = m[1, 0]
        exact_extract[2*i] = exact[ti, si]
        exact_extract[2*i+1] = exact[si, ti]
        ted_extract[2*i] = ted[ti, si]
        ted_extract[2*i+1] = ted[si, ti]
    np.savetxt(dirname + "/MFPT_supernode_estimate.dat", estimate)
    np.savetxt(dirname + "/TED_extract_for_supernode_estimate.dat", ted_extract)
    np.savetxt(dirname + "/MFPT_exact_for_supernode_estimate.dat", exact_extract)

def normalise_by_row(d):
    """Normalise an array row-by-row, ie make each row sum to 1. This
    is useful when creating transition matrices on sub-graphs, or
    randomly-generating transition matrices, or for matrices created
    using a hill-climbing adjustment to transition probabilies: in
    such cases, the transition probabilities are adjusted to take
    account of fitness, so "bad" transitions (to worse fitness) are
    made less likely to be accepted, meaning that each row no longer
    sums to 1. By renormalising we get the true probability after (if
    necessary) multiple rounds of rejection and finally one
    acceptance."""
    for i in range(len(d)):
        sumval = np.sum(d[i])
        d[i] *= 1.0 / sumval
    check_row_sums(d)
    return d

def make_random_matrix(n):
    """Make a random transition matrix on n points."""
    tm = np.random.random((n, n))
    return normalise_by_row(tm)

def make_absorbing(tm, dest):
    """Given a transition matrix, disallow transitions away from the
    given destination -- ie make it an absorbing matrix."""
    for j in range(len(tm)):
        tm[dest, j] = 0.0
    tm[dest, dest] = 1.0
    return tm

def read_transition_matrix(filename):
    """Read a transition matrix from a file and return. The matrix
    will have been written in the right format by some Java code."""
    d = np.genfromtxt(filename)
    return d

def check_row_sums(d):
    # check that each row sums to 1, since each row is the
    # out-probabilities from a single individual. Allow a small margin
    # of error.
    epsilon = 0.000001
    for i, v in enumerate(d):
        if not (1.0 - epsilon < sum(v) < 1.0 + epsilon):
            print("Row doesn't sum to 1: " + str(i) + " " + str(sum(v)))
            assert False
    

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
    # FIXME this is the right method -- avoids python loop
    # x -= np.diag(x) * np.eye(len(x))
        

def get_mfpt(x):
    """Calculate mean-first-passage time of a given transition
    matrix. Set self-transitions to zero. Note that the pysal code
    (ergodic.py) calls it "first-mean-passage-time"."""
    # NB! The ergodic code expects a matrix, not a numpy array. Breaks
    # otherwise.
    x = np.matrix(x)
    x = np.array(ergodic.fmpt(x))
    set_self_transition_zero(x)
    return x
    
def test_matrix_size(n):
    """Test how big the tm can be before get_mfpt becomes
    slow. n = 4000 is fine, n = 10000 starts paging out (at least 30
    minutes)."""
    d = make_random_matrix(n)
    mfpt = get_mfpt(d)
    print("min", np.min(mfpt))
    print("max", np.max(mfpt))

def invert_probabilities(adj):
    """Convert a probability into an edge traversal cost. It's ok to
    ignore runtime floating point problems here, because they're only
    due to zero-probabilities, which result in infinite edge-traversal
    costs, which is what we want. Restore the "raise" behaviour
    afterward."""
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
    traversal costs. See invert_ and deinvert_probabilities."""
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
    paging out (at least 10 minutes) on my machine."""
    # from
    # http://www.depthfirstsearch.net/blog/2009/12/03/computing-all-shortest-paths-in-python/
    n = len(adj)
    for k in range(n):
        adj = np.minimum(adj, np.add.outer(adj[:,k],adj[k,:]))
    return adj


def floyd_warshall_nsteps(adj):
    """Disregard the transition probabilities, other than to see
    whether an edge traversal is allowed or not. Calculate the number
    of steps required to get from each point to each other."""
    x = discretize_probabilities(adj)
    x = floyd_warshall(x)
    set_self_transition_zero(x)
    return x
    
def discretize_probabilities(d):
    """Set the edge cost to 1 if there is a nonzero probability, and
    to infinity if there is a zero probability."""
    retval = np.ones_like(d, dtype=float)
    inf = np.infty
    for i in range(len(d)):
        for j in range(len(d)):
            if not d[i, j] > 0.0:
                retval[i, j] = inf
    return retval

def get_dtp(t):
    """Get D_TP, the distance based on transition probability. D_TP(x,
    y) is the log of TP(x, y), for x != y."""
    x = invert_probabilities(t)
    set_self_transition_zero(x)
    return x

def get_symmetric_version(m):
    """Given an asymmetric matrix, return the symmetrized version,
    which is the mean of the matrix and its transpose."""
    return 0.5 * (m + m.T)

def write_symmetric_remoteness(dirname):
    """Read in the D_TP matrix and the MFPT one, and write out the
    symmetric versions."""
    dtp = np.genfromtxt(dirname + "/D_TP.dat")
    mfpt = np.genfromtxt(dirname + "/MFPT.dat")
    sdtp = get_symmetric_version(dtp)
    ct = get_symmetric_version(mfpt)
    # SD_TP is symmetric transition probability distance
    np.savetxt(dirname + "/SD_TP.dat", sdtp)
    # CT stands for commute time
    np.savetxt(dirname + "/CT.dat", ct)

def get_steady_state(tp):
    """Given a transition probability matrix, use ergodic.steady_state
    to calculate the long-run steady-state, which is a vector
    representing how long the system will spend in each state in the
    long run. If not uniform, that is a bias imposed by the operator
    on the system."""
    import ergodic
    ss = np.array(ergodic.steady_state(np.matrix(tp)))
    ss = np.real(ss) # discard zero imaginary parts
    ss = ss.T[0] # not sure why it ends up with an extra unused dimension
    return ss

def get_Boley_undirected(tp):
    """Boley et al define an undirected graph which "corresponds to" a
    directed graph. Its adjacency matrix is G**s = (Pi * P + P' *
    Pi)/2, where Pi is the steady-state set out along a diagonal and P
    is the transition probability matrix. But we will get the
    transition probability matrix:

    P**s = (P + inv(Pi) * P.T * Pi) / 2

    Note that this matrix is not necessarily symmetric.
    """
    from numpy import dot as d
    P = tp
    Pi = np.diag(get_steady_state(tp))
    return (P + d(d(np.linalg.inv(Pi), P.T), Pi)) / 2.0

def is_symmetric(x):
    return np.allclose(x, x.T)

def Von_Luxburg_amplified_commute_wrapper(dirname):
    # The following requires a symmetric adjacency matrix so we read
    # TP.dat and symmetrize it. can't use the CT.dat which has been
    # written-out, because that corresponds to TP.dat. have to
    # re-calculate CT (done inside the function).
    tp = np.genfromtxt(dirname + "/TP.dat")
    stp = (tp + tp.T) / 2.0
    ct_amp = Von_Luxburg_amplified_commute(stp)
    outfilename = dirname + "/CT_amp.dat"
    np.savetxt(outfilename, ct_amp)

def get_commute_distance_using_Laplacian(S):
    assert np.allclose(S, S.T)

    n = S.shape[0]
    # get laplacian L
    d = np.sum(S,1)
    D = np.diag(d)
    L = (D - S)
    dinv = 1./ np.sum(S, 0)
    
    Linv = linalg.inv(L + np.ones(L.shape)/n) - np.ones(L.shape)/n

    Linv_diag = np.diag(Linv).reshape((n, 1))
    Rexact = Linv_diag * np.ones((1, n)) + np.ones((n, 1)) * Linv_diag.T - 2 * Linv
    
    # convert from a resistance distance to a commute time distance
    vol = np.sum(S)
    Rexact *= vol
    
    return Rexact

def Von_Luxburg_amplified_commute(A):
    R = get_commute_distance_using_Laplacian(A)

    n = A.shape[0]
    
    # compute commute time limit expression:
    d = np.sum(A, 1)    
    Rlimit = np.tile((1. / d), (n, 1)).T + np.tile((1. / d), (n, 1))

    # compute correction term u_{ij}: 
    tmp = np.tile(np.diag(A), (n, 1)).T / np.tile(d, (n, 1))
    tmp2 = np.tile(d, (n, 1)).T * np.tile(d, (n, 1))
    uAu = tmp + tmp.T - 2 * A / tmp2

    # compute amplified commute: 
    tmp3 = R - Rlimit - uAu

    # enforce 0 diagonal: 
    D = tmp3 - np.diag(np.diag(tmp3))

    return D

def read_and_get_Von_Luxburg_approximations(dirname):
    # assumes TP and MFPT have been calculated and written out already
    t = read_transition_matrix(dirname + "/TP.dat")
    mfpt = np.genfromtxt(dirname + "/MFPT.dat")
    ones = np.ones_like(t)
    n = t.shape[0]
    
    # d_v is the in-degree, which is the column sum. this
    # multiplication does the right thing
    d_v = np.sum(t, axis=0) * ones
    # transpose to get d_u
    d_u = d_v.T
    # vol(G) is sum of degrees, which we interpret in the directed
    # case as sum of out-degrees (which is anyway equal to sum of
    # in-degrees); because weights are transition probabilities, each
    # out-degree = 1; so vol(G) = n.
    vol_G = float(n)
    
    mfpt_vla = vol_G * 1.0 / d_v
    set_self_transition_zero(mfpt_vla)
    outfilename = dirname + "/MFPT_VLA.dat"
    np.savetxt(outfilename, mfpt_vla)

    ct_vla = vol_G * (1.0 / d_u + 1.0 / d_v)
    set_self_transition_zero(ct_vla)
    outfilename = dirname + "/CT_VLA.dat"
    np.savetxt(outfilename, ct_vla)
    

def Laplacian_matrix(A):
    """Copied from
    http://networkx.lanl.gov/_modules/networkx/linalg/laplacianmatrix.html. Eg

>>> A = np.array([[0, 1, 0, 0, 1, 0],
                  [1, 0, 1, 0, 1, 0],
                  [0, 1, 0, 1, 0, 0],
                  [0, 0, 1, 0, 1, 1], 
                  [1, 1, 0, 1, 0, 0],
                  [0, 0, 0, 1, 0, 0]])
>>> Laplacian_matrix(A)
array([[ 2., -1.,  0.,  0., -1.,  0.],
       [-1.,  3., -1.,  0., -1.,  0.],
       [ 0., -1.,  2., -1.,  0.,  0.],
       [ 0.,  0., -1.,  3., -1., -1.],
       [-1., -1.,  0., -1.,  3.,  0.],
       [ 0.,  0.,  0., -1.,  0.,  1.]])
"""
    I=np.identity(A.shape[0])
    D=I*np.sum(A,axis=1)
    L=D-A
    return L

def Yen_correction(m):
    """From Yen et al,
    http://link.springer.com/content/pdf/10.1007%2F978-3-540-71701-0_117.pdf,
    referred to by Von Luxburg et al."""
    L = np.linalg.pinv(Laplacian_matrix(m))
    sigma = np.std(L) # a normalising factor
    a = 7.0 # The value used by Yen et al, determined experimentally
    K_Yen = 1.0 / (1 + a * np.exp(-L / sigma))
    return K_Yen

def Brand_correction(m):
    """Placeholder. From Von Luxburg et al."""
    # assumes mfpt has been calculated and written out already
    # mfpt = np.genfromtxt(dirname + "/MFPT.dat")
    # R = mfpt / vol_G
    n = m.shape[0]
    sum_all = np.sum(R)
    sum_ik_kj = np.zeros_like(R)

    # this is far too slow -- need to reformulate in matrix terms
    for i in range(n):
        for j in range(n):
            sum_ik_kj = sum(R[i][k] + R[k][j] for k in range(n))
    K = 1/2.0 * (-R + 1./n * (sum_ik + sum_kj) - 1.0/(n*n) * sum_all)
    K_B = K / np.sqrt(Kii * Kjj) # (need to express KiiKjj as a matrix)
    return K_B
    
def read_and_get_dtp_mfpt_sp_steps(dirname):

    if os.path.exists(dirname + "/TP_nonnormalised.dat"):
        t = read_transition_matrix(dirname + "/TP_nonnormalised.dat")
        # it's a non-normalised matrix, possibly representing a
        # uniformly sampled sub-graph, a hill-climb-sampled graph, or
        # similar.
        t = normalise_by_row(t)
        outfilename = dirname + "/TP.dat"
        np.savetxt(outfilename, t)
    else:
        t = read_transition_matrix(dirname + "/TP.dat")
        check_row_sums(t)

    # This gets D_TP, which is just the transition probability inverted
    d = get_dtp(t)
    outfilename = dirname + "/D_TP.dat"
    np.savetxt(outfilename, d)
    
    # This gets the mean first passage time, ie the expected length of
    # a random walk.
    f = get_mfpt(t)
    outfilename = dirname + "/MFPT.dat"
    np.savetxt(outfilename, f)

    # This gets the cost of the shortest path between pairs. The cost
    # of an edge is the negative log of its probability.
    h = floyd_warshall_probabilities(t)
    outfilename = dirname + "/SP.dat"
    np.savetxt(outfilename, h)

    # this gets the minimum number of steps required to go between
    # pairs, disregarding probabilities. Only interesting if some
    # edges are absent (ie edge probability is zero).
    p = floyd_warshall_nsteps(t)
    outfilename = dirname + "/STEPS.dat"
    np.savetxt(outfilename, p)

def hamming_distance(x, y):
    return np.sum(x != y)
    
def generate_ga_tm(dirname, pmut=None):
    """For a bitstring (genetic algorithm) representation of a given
    length, generate a transition matrix with the mutation probability
    pmut. Also generate the Hamming distances. If pmut=None (default),
    exactly one bitflip is performed per individual, rather than using
    a per-gene mutation probability."""
    
    length = int(dirname.strip("/").split("_")[2])
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
    outfilename = dirname + "/TP.dat"
    np.savetxt(outfilename, tm)
    outfilename = dirname + "/Hamming.dat"
    np.savetxt(outfilename, hm)

def land_of_oz_matrix():
    """From Kemeny & Snell 1976. The states are rain, nice and
    snow."""
    return np.array([[.5, .25, .25], [.5, .0, .5], [.25, .25, .5]])

def simulate_random_walk(f, nsteps, selected, nsaves):
    """f is the transition function. nsteps is the number of steps.
    selected is the set of states (can be a list of integers) in which
    we're interested. nsaves is the max number of samples to save for
    each pair. Return an array of samples for MFPT for each pair (i,
    j) from selected."""

    n = len(selected)
    # samples is all nans to start
    samples = np.nan + np.zeros((n, n, nsaves))
    # in the rw_started matrix, the (i,j)th entry is the time-step at
    # which a walk from i to j began. -1 means we are not currently in
    # a walk from i to j. can't be in a walk from i to j and from j to
    # i simultaneously.
    rw_started = -1 * np.ones((n, n), dtype='int64')
    # for each transition there are nsaves spaces to save samples --
    # this says where to save
    saveidx = np.zeros((n, n), dtype='int')

    def print_state(v, rw_started, samples, t):
        print("%d: %d" % (t, v))

    # su and sv are states (may be integers). u and v are the indices
    # into selected of su and sv
    sv = selected[0]
    for t in range(nsteps):
        print_state(sv, rw_started, samples, t)
        
        # when we see a state sv which is of interest
        if sv in selected:
            v = selected.index(sv)
            
            for u, su in enumerate(selected):

                if rw_started[u,v] == rw_started[v,u] == -1:
                    rw_started[v,u] = t

                elif rw_started[u,v] > -1:
                    # we are currently in a walk from u to v, and we
                    # have just reached v, so save a sample (now -
                    # start of walk), if there is space
                    if saveidx[u,v] < nsaves:
                        samples[u,v,saveidx[u,v]] = (t - rw_started[u,v])
                        # update the space left
                        saveidx[u,v] += 1
                    # now we start a walk from v to u. because of the
                    # order of these two assignments, this does the
                    # right thing if v == u.
                    rw_started[u,v] = -1
                    rw_started[v,u] = t

                elif rw_started[v,u] > -1:
                    # we are in a walk from v to u, and have
                    # re-encountered v. ignore.
                    pass

                else:
                    raise InternalError
                    
        # transition to a new state
        sv = f(sv)
    return samples

def roulette_wheel(a):
    """Randomly sample an index, weighted by the given sequence of
    probabilities."""
    s = np.sum(a)
    assert s > 0
    r = random.random() * s
    cumsum = 0.0
    for i, v in enumerate(a):
        cumsum += v
        if r < cumsum:
            return i
    raise ValueError("Unexpected: failed to find a slot in roulette wheel")

    
def test_random_walk():
    oz = land_of_oz_matrix()
    samples = simulate_random_walk(
        lambda i: roulette_wheel(tp[i]), # transition according to tp
        200, [0, 1, 2], 10)
    print(samples)
    mfpte = scipy.stats.nanmean(samples, axis=2)
    mfpte_std = scipy.stats.nanstd(samples, axis=2)
    mfpt = get_mfpt(oz)
    print("MFPTE")
    print(mfpte)
    print("MFPTE_STD")
    print(mfpte_std)
    print("MFPT")
    print(mfpt)
    
if __name__ == "__main__":
    dirname = sys.argv[1]

    # if "depth" in dirname:
    #     # Matrices have already been generated by Java code.
    #     pass
    # elif "per_ind" in dirname:
    #     generate_ga_tm(dirname)
    # else:
    #     generate_ga_tm(dirname, 0.1)
    # read_and_get_dtp_mfpt_sp_steps(dirname)
    # write_symmetric_remoteness(dirname)
    # read_and_get_Von_Luxburg_approximations(dirname)
    # estimate_MFPT_with_supernode(dirname)
    # analyse_random_walk(dirname)
    test_random_walk()
