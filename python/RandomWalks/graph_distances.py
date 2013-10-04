#!/usr/bin/env python

"""This module provides methods related to various novel graph
distances, from research by Boley etal, Von Luxburg etal, Kivimaki
etal."""

import numpy as np
import scipy.stats
import scipy.linalg as linalg
from numpy import dot as d
import random
import sys

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
    """Original code copyright (C) Ulrike Von Luxburg, Python
    implementation by me (James McDermott)."""

    n = S.shape[0]
    L = Laplacian_matrix(S)
    dinv = 1./ np.sum(S, 0)
    
    Linv = linalg.inv(L + np.ones(L.shape)/n) - np.ones(L.shape)/n

    Linv_diag = np.diag(Linv).reshape((n, 1))
    Rexact = Linv_diag * np.ones((1, n)) + np.ones((n, 1)) * Linv_diag.T - 2 * Linv
    
    # convert from a resistance distance to a commute time distance
    vol = np.sum(S)
    Rexact *= vol
    
    return Rexact

def Von_Luxburg_amplified_commute(A):
    """From Von Luxburg etal, "Getting lost in space: Large sample
    analysis of the commute distance". Original code copyright (C)
    Ulrike Von Luxburg, Python implementation by me (James
    McDermott)."""
    
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
    """From Von Luxburg etal, 2011, "Hitting and commute times in
    large graphs are often misleading"."""
    
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
    assert is_symmetric(A)
    I=np.identity(A.shape[0])
    D=I*np.sum(A,axis=1)
    L=D-A
    return L
    

def RSP_and_FE_distances(A, beta, C=None):
    """Calculate the randomised shortest path distance and free-energy
    distance, as defined in "Developments in the theory of randomized
    shortest paths with a comparison of graph node distances",
    Kivim\"{a}ki, Shimbo, and Saerens, Physica A, 2013.

    Arguments
    
    - A:     Adjacency matrix, whose elements represent affinities between
             nodes, which define the reference transition
             probabilities. A can be asymmetric (for directed graphs).
             Distances between nodes that are not strongly connected
             are Inf.
             
    - beta   beta should lie more or less between 10^-8 and 20, but this 
             depends on the size of the graph and the magnitude of the
             costs. When beta --> 0, we obtain the commute cost distances.
             When beta --> \infty, we obtain the shortest path (lowest
             cost) distances.
             
    - C      Cost matrix, whose elements represent the cost of traversing
             an edge of the graph. Infinite costs can be marked as zeros
             (zero costs are anyway not allowed). If C is not provided,
             then the costs will be set by default as c_ij = 1/a_ij.
    
    Returns D_RSP: the RSP dissimilarity matrix
            D_FE:  the free energy distance matrix
    
    Original Matlab code and comments (c) Ilkka Kivim\"{a}ki 2013
    
    Transliterated to Python/Numpy by James McDermott
    <jamesmichaelmcdermott@gmail.com>. Helpful guides to this type of
    transliteration:
    http://mathesaurus.sourceforge.net/matlab-numpy.html,
    http://wiki.scipy.org/NumPy_for_Matlab_Users,
    http://wiki.scipy.org/Tentative_NumPy_Tutorial
    """

    max = np.finfo('d').max
    eps = 0.00000001

    # If A is integer-valued, and beta is floating-point, can get an
    # error in the matrix inversion, so convert A to float here. I
    # can't explain why beta being floating-point is related to the
    # problem. Anyway, this also converts in case it was a matrix, or
    # was sparse.
    A = np.array(A, dtype=np.float)

    A[A < eps] = 0.0
    n, m = A.shape
    if n != m:
        raise ValueError("The input matrix A must be square")

    if C is None:
        C = A.copy()
        C[A >= eps] = 1.0/A[A >= eps]
        C[A < eps] = max

    # check beta value?
    if beta < eps or beta > 20.0:
        raise ValueError("The value for beta is outside the expected range, 0 to 20.0")

    ones = np.ones(n)
    onesT = np.ones((n, 1))
    I = np.eye(n)

    # Computation of Pref, the reference transition probability matrix
    tmp = A.copy()
    s = np.sum(tmp, 1)
    s[s == 0] = 1 # avoid zero-division
    Pref = tmp / (s * onesT).T

    # Computation of the W and Z matrices
    W = np.exp(-beta * C) * Pref

    # compute Z
    Z = linalg.inv(I - W)

    # Computation of Z*(C.*W)*Z avoiding zero-division errors:
    numerator = d(d(Z, (C * W)), Z)
    D_nonabs = np.zeros((n, n))
    
    indx = (numerator > 0) & (Z > 0)
    D_nonabs[indx] = numerator[indx] / Z[indx]
    D_nonabs[~indx] = np.infty
    # D_nonabs above actually gives the expected costs of non-hitting paths
    # from i to j.

    # Expected costs of hitting paths -- avoid a possible inf-inf
    # which can arise with isolated nodes and would give a NaN -- we
    # prefer to have inf in that case.
    C_RSP = np.zeros((n, n))
    diag_D = d(onesT, np.diag(D_nonabs).reshape((1, n)))
    indx = ~np.isinf(diag_D)
    C_RSP[indx] = D_nonabs[indx] - diag_D[indx]
    C_RSP[~indx] = np.infty

    # symmetrization
    D_RSP = 0.5 * (C_RSP + C_RSP.T)

    # Free energies and symmetrization:
    Dh_1 = np.diag(1.0/np.diag(Z))
    Zh = d(Z, Dh_1)

    # If there any 0 values in Zh (because of isolated nodes), taking
    # log will raise a divide-by-zero error -- ignore it
    np.seterr(divide='ignore')
    FE = -np.log(Zh)/beta
    np.seterr(divide='raise')
    D_FE = 0.5 * (FE + FE.T)

    # Just in case, set diagonals to zero:
    np.fill_diagonal(D_RSP, 0.0)
    np.fill_diagonal(D_FE, 0.0)

    return D_RSP, D_FE

def test_Kivimaki():
    """From Kivimaki code (note beta = 1.0):
octave-3.4.0:11> A = [0 0 0 0; 0 0 0 1; 0 0 0 1; 0 1 1 0]
A =

   0   0   0   0
   0   0   0   1
   0   0   0   1
   0   1   1   0

octave-3.4.0:12> D_RSP, D_FE = RSP_distances(A, 1)
error: `D_RSP' undefined near line 12 column 1
octave-3.4.0:12> [D_RSP, D_FE] = RSP_distances(A, 1)
warning: The distances have Inf values! This is either because the graph is not connected or beta is too small or too large
D_RSP =

   0.00000       Inf       Inf       Inf
       Inf   0.00000   2.14516   1.07258
       Inf   2.14516   0.00000   1.07258
       Inf   1.07258   1.07258   0.00000

D_FE =

   0.00000       Inf       Inf       Inf
       Inf   0.00000   2.62308   1.31154
       Inf   2.62308   0.00000   1.31154
       Inf   1.31154   1.31154   0.00000

octave-3.4.0:273> A = [.5 .25 .25; .25 .0 .75; .25 .75 .0]
A =

   0.50000   0.25000   0.25000
   0.25000   0.00000   0.75000
   0.25000   0.75000   0.00000

octave-3.4.0:274> [D_RSP, D_FE] = RSP_distances(A, 1.)
D_RSP =

   0.00000   4.34699   4.34699
   4.34699   0.00000   1.33429
   4.34699   1.33429   0.00000

D_FE =

   0.00000   5.15091   5.15091
   5.15091   0.00000   1.62088
   5.15091   1.62088   0.00000

   """       
    A = np.array([[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 1], [0, 1, 1, 0]])
    beta = 1.0
    D_RSP, D_FE = RSP_and_FE_distances(A, beta)
    print(D_RSP)
    print(D_FE)

    A = np.array([[.5, .25, .25], [.25, .0, .75], [.25, .75, .0]])
    beta = 1.0
    D_RSP, D_FE = RSP_and_FE_distances(A, beta)
    print(D_RSP)
    print(D_FE)

    # A = np.array([[0.5, 0.5], [0.1, 0.9]])


def get_commute_distance_using_Laplacian(S):
    assert np.allclose(S, S.T)

    n = S.shape[0]
    # compute global graph laplacian: 
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
    # FIXME should be able to use our existing MFPT and scale by the
    # right constant to get R:
    
    # mfpt = get_mfpt(A)
    # R = mfpt + mfpt.T
    # # convert our commute time to a resistance distance as used by VL
    # R /= vol 

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

def test_Von_Luxburg():
    vl_matlab_result = """
    octave-3.4.0:201> b
b =

   0.50000   0.25000   0.25000
   0.25000   0.00000   0.75000
   0.25000   0.75000   0.00000

octave-3.4.0:202> compute_amplified_commute_distance(b)
ans =

   0.00000   4.85714   4.85714
   4.85714   0.00000   2.92857
   4.85714   2.92857   0.00000
"""
    b = np.array([[.5, .25, .25], [.25, .0, .75], [.25, .75, .0]])
    print("b")
    print(b)
    print("Laplacian:")
    print(get_commute_distance_using_Laplacian(b))
    print("Result from VL Matlab code:")
    print(vl_matlab_result)
    print("Our result:")
    print Von_Luxburg_amplified_commute(b)


if __name__ == "__main__":
    test_Kivimaki()
    test_Von_Luxburg()
