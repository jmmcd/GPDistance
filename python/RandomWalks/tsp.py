#!/usr/bin/env python

import random
import collections
import itertools

###################################################################
# TSP stuff
###################################################################

def swap_two(p, ab=None):
    """Swap-two just swaps any two elements of a permutation. As long as
    they are distinct a new permutation is formed. We canonicalise
    after. Can pass in ab in which case the mutation is deterministic,
    useful for analysis."""
    if ab is None:
        a, b = sorted(random.sample(xrange(len(p)), 2))
    else:
        a, b = ab
    p = p[:]
    p[a], p[b] = p[b], p[a]
    return canonicalise(p)

def swap_adj(p, ab=None):
    """Choose an element and swap with the following element (if we choose
    the last, swap with the first). Canonicalise after. Can pass in ab
    in which case the mutation is deterministic, useful for
    analysis."""
    if ab is None:
        n = len(p)
        a = random.randrange(n)
        b = (a+1) % n
    else:
        a, b = ab
    p = p[:]
    p[a], p[b] = p[b], p[a]
    return canonicalise(p)

def two_opt(p):
    """2-opt means choosing any two non-contiguous edges ab and cd,
    chopping them, and then reconnecting (such that the result is
    still a complete tour). There are actually two ways of doing it --
    one is the identity, and one gives a new tour."""
    n = len(p)
    a = random.randrange(n)
    b = random.randrange(n-2) # b=a-1, b=a, and b=a+1 are invalid
    if (b % n) in ((a-1) % n, a % n, (a+1) % n):
        b = (b+3) % n
    a, b = min(a, b), max(a, b)
    p = p[:a+1] + p[b:a:-1] + p[b+1:]
    return canonicalise(p)

def three_opt_broad(p):
    return three_opt(p, broad=True)

def three_opt(p, broad=False):
    """In the broad sense, 3-opt means choosing any three edges ab, cd
    and ef and chopping them, and then reconnecting (such that the
    result is still a complete tour). There are eight ways of doing
    it. One is the identity, 3 are 2-opt moves (because either ab, cd,
    or ef is reconnected), and 4 are 3-opt moves (in the narrower
    sense)."""
    n = len(p)
    # choose 3 unique edges defined by their first node
    a, c, e = random.sample(range(n+1), 3)
    # without loss of generality, sort
    a, c, e = sorted([a, c, e])
    b, d, f = a+1, c+1, e+1

    if broad == True:
        which = random.randint(0, 7) # allow any of the 8
    else:
        which = random.choice([3, 4, 5, 6]) # allow only strict 3-opt

    # in the following slices, the nodes abcdef are referred to by
    # name. x:y:-1 means step backwards. anything like c+1 or d-1
    # refers to c or d, but to include the item itself, we use the +1
    # or -1 in the slice
    if which == 0:
        sol = p[:a+1] + p[b:c+1]    + p[d:e+1]    + p[f:] # identity
    elif which == 1:
        sol = p[:a+1] + p[b:c+1]    + p[e:d-1:-1] + p[f:] # 2-opt
    elif which == 2:
        sol = p[:a+1] + p[c:b-1:-1] + p[d:e+1]    + p[f:] # 2-opt
    elif which == 3:
        sol = p[:a+1] + p[c:b-1:-1] + p[e:d-1:-1] + p[f:] # 3-opt
    elif which == 4:
        sol = p[:a+1] + p[d:e+1]    + p[b:c+1]    + p[f:] # 3-opt
    elif which == 5:
        sol = p[:a+1] + p[d:e+1]    + p[c:b-1:-1] + p[f:] # 3-opt
    elif which == 6:
        sol = p[:a+1] + p[e:d-1:-1] + p[b:c+1]    + p[f:] # 3-opt
    elif which == 7:
        sol = p[:a+1] + p[e:d-1:-1] + p[c:b-1:-1] + p[f:] # 2-opt

    return canonicalise(sol)

def canonicalise(p):
    """In any TSP, 01234 is equivalent to 23410. We canonicalise on
    the former. In a symmetric TSP, 01234 is equivalent to 04321. We
    canonicalise on the former. The general recipe is: p[0] is 0, and
    p[1] is the smaller of p[1]'s neighbours (so p[-1] is the
    larger)."""
    if p[0] != 0:
        i = p.index(0)
        p[:] = p[i:] + p[:i]
    if p[1] > p[-1]:
        p[1:] = p[-1:0:-1]
    return p


def tsp_tours(n):
    """Generate all tours of length n. A tour is a permutation. But we
    canonicalise as above."""
    for p in itertools.permutations(range(n)):
        if p[0] != 0: continue
        if p[1] > p[-1]: continue
        yield p

def sample_transitions(n, move="2opt", nsamples=10000):
    length = len(list(tsp_tours(n)))
    tm = np.zeros((length, length))
    delta = 1.0 / nsamples

    if move == "3opt":
        move = three_opt
    elif move == "3opt_broad":
        move = three_opt_broad
    elif move == "2opt":
        move = two_opt
    elif move == "swap":
        move = swap_two
    elif move == "swap_adj":
        move = swap_adj
    else:
        raise ValueError("Unexpected move type: " + str(move))

    tours_to_ints = {}
    for i, tour in enumerate(tsp_tours(n)):
        tours_to_ints[tour] = i

    for i, tour in enumerate(tsp_tours(n)):
        for j in range(nsamples):
            t = move(list(tour))
            tm[i][tours_to_ints[tuple(t)]] += delta
    return tm

def kendall_tau_permutation_distance(t1, t2):
    corr, p = scipy.stats.kendalltau(t1, t2)
    return 1.0 - corr

def kendall_tau_permutation_distances(n):
    m = len(list(tsp_tours(n)))
    kt = np.zeros((m, m))
    for i, ti in enumerate(tsp_tours(n)):
        for j, tj in enumerate(tsp_tours(n)):
            kt[i][j] = kendall_tau_permutation_distance(ti, tj)
    return kt

def tsp_tm_wrapper(dirname, move="2opt"):
    # dirname should be <dir>/tsp_length_6_2opt, for example
    t = dirname.find("tsp_length_")
    length = dirname[t:].split("_")[2]
    length = int(length)
    tm = sample_transitions(length, move)
    outfilename = dirname + "/TP.dat"
    np.savetxt(outfilename, tm)
    kt = kendall_tau_permutation_distances(length)
    outfilename = dirname + "/KendallTau.dat"
    np.savetxt(outfilename, kt)

#
# end of TSP stuff
#######################################################


def test_2opt():
    n = 10
    m = 100000
    x = y = z = 0
    for i in range(m):
        a = random.randrange(n)
        
        # r = range(n)
        # r.remove(a)
        # r.remove((a+1) % n)
        # r.remove((a-1) % n)
        # print r
        # b = random.choice(r)

        # if a == 0:
        #     x += 1
        #     b = random.randrange(2, n-1)
        # elif a == n-1:
        #     y += 1
        #     b = random.randrange(1, n-2)
        # else:
        #     z += 1
        #     b = random.choice(range(0, a-1) + range(a+2, n))
            
        b = random.randrange(n-3)
        print a, b
        if b in ((a-1) % n, a, (a+1) % n):
            b = (b+3) % n
        print a, b
        a, b = sorted([a, b])
        print a, b
        print
        yield (a, b)
    print float(x) / m, float(y) / m, float(z) / m
    c = n * (n-3) / 2
    print c
    print m / float(c)

print collections.Counter(test_2opt()).most_common(35)
