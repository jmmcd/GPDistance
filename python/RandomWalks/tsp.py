#!/usr/bin/env python

import random
import collections
import itertools

###################################################################
# TSP stuff
#
# for general TSP (no assumption of symmetric costs)
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

def two_opt(p, ab=None):
    """2-opt means choosing any two non-contiguous edges ab and cd,
    chopping them, and then reconnecting (such that the result is
    still a complete tour). There are actually two ways of doing it --
    one is the identity, and one gives a new tour.

    There are n ways of choosing a, then n-3 ways of choosing c, but
    because order is unimportant between these the number of ways to
    choose a and c is n(n-3)/2
    """
    if ab is None:
        n = len(p)
        a = random.randrange(n)
        # b=a-1, b=a, and b=a+1 are invalid
        b = random.randrange(a+2, a+n-1) % n
        a, b = min(a, b), max(a, b)
    else:
        a, b = ab
    p = p[:a+1] + p[b:a:-1] + p[b+1:]
    return canonicalise(p)

def twoh_opt(p, abc=None):
    """2h-opt aka 2.5-opt. Choose AB and C, and move node C to form ACB.
    Bentley: "Fast Algorithms for Geometric Traveling Salesman Problems"

    There are n ways of choosing a and n-2 ways of choosing c, so
    n(n-2) altogether. But some of these turn out to be equivalent, eg
    (0, 1, 2, 3, 4). Choose A=2, C=1, get (0, 2, 1, 3, 4). Choose A=0,
    C=2, again get (0, 2, 1, 3, 4). But there's only one way to get
    (0, 2, 3, 1, 4). Hence some neighbours are more likely than
    others.

    """
    if abc is None:
        n = len(p)
        a = random.randrange(n)
        b = (a+1) % n
        c = random.randrange(a+2, a+n) % n # c can be any node other than a or b
    else:
        a, b, c = abc
    p = p[:]
    p.insert(a+1, p[c])
    if a < c:
        del p[c+1]
    else:
        del p[c]
    return canonicalise(p)

def three_opt_broad(p):
    return three_opt(p, broad=True)

def three_opt(p, broad=False):
    """In the broad sense, 3-opt means choosing any three non-contiguous
    edges ab, cd and ef and chopping them, and then reconnecting (such
    that the result is still a complete tour).

    There are n choices for the first edge.

    For the second edge, choices which touch the first edge are
    disallowed. There are 2 ways of choosing the second edge such that
    it's within 2 of the first edge, and thereafter there are n-5 ways
    for the third.
    
    There are n-5 ways of choosing the second edge such that it's
    *not* within 2 of the first edge, and thereafter n-6 ways for the
    third.
    
    Hence the total number of choices of ab cd ef is n*2*(n-5) +
    n*(n-5)*(n-6). But the order of a, c, e is unimportant, so we
    divide by 3! = 6, so the expression is n(n-4)(n-5)/6.

    Eg with a 6-node tour there are just two ways to choose the three
    edges -- either start on 0 and take every second one, or start on
    1 and take every second one. We use a special case for this
    because otherwise the wrong choice of c makes no choice of e
    valid. 

    Given the choice of ab cd ef, there are eight ways of
    reconnecting. One is the identity, 3 are 2-opt moves (because
    either ab, cd, or ef is reconnected), and 4 are 3-opt moves (in
    the narrower sense).

    Therefore I think the total number of options, all equally
    probable, is 2n(n-4)(n-5)/3!

    """
    n = len(p)
    
    # choose 3 unique non-contiguous edges defined by their first node
    if n <= 5:
        raise ValueError
    elif n == 6:
        a = random.randrange(2)
        b, c, d, e, f = range(a+1, a+6)
    else:

        # FIXME this still isn't right
        a = random.randrange(n)
        c = random.randrange(a+2, a+n-1) % n
        print "a, c", a, c
        r = range(n)
        r.remove((a-1) % n)
        r.remove(a)
        r.remove((a+1) % n)
        print (c-1)%n, r
        try:
            r.remove((c-1) % n)
        except:
            print "failed to remove", (c-1)%n
        print c, r
        r.remove(c)
        print (c+1)%n, r
        try:
            r.remove((c+1) % n)
        except:
            print "failed to remove", (c+1)%n
        print r
        e = random.choice(r)

        a, c, e = sorted([a, c, e])
        b = (a+1) % n
        d = (c+1) % n
        f = (e+1) % n

    print "abcdef", a, b, c, d, e, f

    if broad == True:
        which = random.choice([1, 2, 3, 4, 5, 6, 7]) # allow any of the 2-opt or 3-opt
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
    """In any TSP, 01234 is equivalent to 23401. We canonicalise on the
    former. In a symmetric TSP, 01234 is equivalent to 04321. This
    code is for a non-symmetric TSP. """
    if p[0] != 0:
        i = p.index(0)
        p[:] = p[i:] + p[:i]
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

def test_op(op):
    n = 7
    m = 100000
    for i in range(m):

        p = range(n)
        p = op(p)
        yield tuple(p)

print collections.Counter(test_op(twoh_opt)).most_common(100)

