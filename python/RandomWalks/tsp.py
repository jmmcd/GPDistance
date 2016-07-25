#!/usr/bin/env python

import random
import collections
import itertools
import numpy as np
import scipy.stats, scipy.misc
import math

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

def _three_opt_choose_edges_unused(n):
    """Choose 3 unique non-contiguous edges defined by their first
    node. This method should not be used, because it does not give a
    uniform sampling of the possible triples.

    """
    if n <= 5:
        raise ValueError
    elif n == 6:
        a = random.randrange(2)
        b, c, d, e, f = a+1, a+2, a+3, a+4, (a+5)%n
    else:

        a = random.randrange(n)
        c = random.randrange(a+2, a+n-1) % n

        a, c = min(a, c), max(a, c) # now a < c, so a < n-2 and c > 2

        minv = 0
        maxv = n-1
        if a == 0:
            maxv = n-2
        if c == n-1:
            minv = 1
        e = random.choice(range(minv, a-1) + range(a+2, c-1) + range(c+2, maxv+1))

        a, c, e = sorted([a, c, e])
        b = a+1
        d = c+1
        f = (e+1)%n

        assert len(set([a, b, c, d, e, f])) == 6

        return a, b, c, d, e, f


def _three_opt_choose_edges_iter(n):
    """Choose 3 unique non-contiguous edges defined by their first node.
    This iterator yields each possible triple in turn."""
    for a in range(n-4):
        for c in range(a+2, n-2):
            if a == 0:
                max_e = n-2
            else:
                max_e = n-1
            for e in range(c+2, max_e+1):
                yield a, a+1, c, c+1, e, (e+1)%n

def _three_opt_choose_edges(n):
    """Choose 3 unique non-contiguous edges defined by their first node.
    This function chooses a random triple, uniformly sampled."""
    a = random.randrange(n-4)
    c = random.randrange(a+2, n-2)
    if a == 0:
        max_e = n-2
    else:
        max_e = n-1
    e = random.randrange(c+2, max_e+1)
    return (a, a+1, c, c+1, e, (e+1)%n)

def _three_opt_iter(n):
    """Choose 3 unique non-contiguous edges defined by their first node.
    This iterator yields each possible triple in turn, but there will
    be duplicates, so it should not be used."""
    for _a in range(n):
        for _c in range(_a+2, n):
            a = _a
            c = _c % n
            a, c = min(a, c), max(a, c) # now a < c, so a < n-2 and c > 2

            minv = 0
            maxv = n-1
            if a == 0:
                maxv = n-2
            if c == n-1:
                minv = 1

            r = range(minv, a-1) + range(a+2, c-1) + range(c+2, maxv+1)
            for e in r:

                b = a+1
                d = c+1
                f = (e+1)%n

                assert len(set([a, b, c, d, e, f])) == 6

                yield sorted([a, b, c, d, e, f])

def three_opt_n_triples(n):
    return n * (n-4) * (n-5) / 6

def three_opt_n_neighbours(n, broad=False):
    # see below for workings
    # note that either (n-4) or (n-5) is divisible by 2
    # and one of n, (n-4) or (n-5) is divisible by 3
    # hence dividing by 6 guarantees an integer result
    if broad:
        return 7 * n * (n-4) * (n-5) / 6
    else:
        return 4 * n * (n-4) * (n-5) / 6

def three_opt_deterministic(p, abcdef, which):
    """Given the edges ab cd and ef to cut, and a choice of which of them
    to reconnect, carry out the cut and reconnection.

    """
    n = len(p)

    a, b, c, d, e, f = abcdef

    # in the following slices, the nodes abcdef are referred to by
    # name. x:y:-1 means step backwards. anything like c+1 or d-1
    # refers to c or d, but to include the item itself, we use the +1
    # or -1 in the slice

    if f == 0:
        minv = 1
        maxv = 0
    else:
        minv = 0
        maxv = n-1
    if which == 0:
        sol = p[minv:a+1] + p[b:c+1]    + p[d:e+1]    + p[f:maxv+1] # identity
    elif which == 1:
        sol = p[minv:a+1] + p[b:c+1]    + p[e:d-1:-1] + p[f:maxv+1] # 2-opt
    elif which == 2:
        sol = p[minv:a+1] + p[c:b-1:-1] + p[d:e+1]    + p[f:maxv+1] # 2-opt
    elif which == 3:
        sol = p[minv:a+1] + p[c:b-1:-1] + p[e:d-1:-1] + p[f:maxv+1] # 3-opt
    elif which == 4:
        sol = p[minv:a+1] + p[d:e+1]    + p[b:c+1]    + p[f:maxv+1] # 3-opt
    elif which == 5:
        sol = p[minv:a+1] + p[d:e+1]    + p[c:b-1:-1] + p[f:maxv+1] # 3-opt
    elif which == 6:
        sol = p[minv:a+1] + p[e:d-1:-1] + p[b:c+1]    + p[f:maxv+1] # 3-opt
    elif which == 7:
        sol = p[minv:a+1] + p[e:d-1:-1] + p[c:b-1:-1] + p[f:maxv+1] # 2-opt

    if len(sol) != n:
        print("bad length:", sol)
        raise
    return canonicalise(sol)

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
    n*(n-5)*(n-6) = n(n-4)(n-5). But the order of a, c, e is
    unimportant, so we divide by 3! = 6, so the expression is
    n(n-4)(n-5)/6. Eg for n=6 we have 2; for n=7 we have 7; for n = 8
    we have 16; for n = 9 we have 30; for n=10, we have 50 choices.

    Eg with a 6-node tour there are just two ways to choose the three
    edges -- either start on 0 and take every second one, or start on
    1 and take every second one. We use a special case for this 6-node
    case because otherwise the wrong choice of c makes no choice of e
    valid.

    Another way to think of it: without loss of generality we can
    assume a < c < e, so we have (n-4) choices for a and then (n-2-a)
    choices for c, and then either (n-c) or (n-c-1) choices for
    e. This is the idea we actually use to generate.

    Another way to think of it: how many "interior" triangles are
    there in the cycle? An interior triangle is one which does not
    share an edge with the cycle. An interior triangle is specified by
    three unique non-adjacent nodes.

    See https://oeis.org/A005581

    Given the choice of ab cd ef, there are eight ways of
    reconnecting. One is the identity, 3 are 2-opt moves (because
    either ab, cd, or ef is reconnected), and 4 are 3-opt moves (in
    the narrower sense).

    Therefore the number of options, all equally probable, is
    7n(n-4)(n-5)/6 (for broad interpretation) or 4n(n-4)(n-5)/6 (for
    strict interpretation).

    """
    n = len(p)

    # choose 3 unique non-contiguous edges
    # this will be inefficient for larger n
    a, b, c, d, e, f = random.choice(list(_three_opt_choose_edges_iter(n)))

    if broad == True:
        # allow any of the 2-opt or 3-opt
        which = random.choice([1, 2, 3, 4, 5, 6, 7])
    else:
        # allow only strict 3-opt
        which = random.choice([3, 4, 5, 6])

    return three_opt_deterministic(p, (a, b, c, d, e, f), which)

def three_opt_broad(p):
    return three_opt(p, broad=True)

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
        yield list(p)

def get_tm(n, move="two_opt"):
    tours = list(tsp_tours(n))
    length = len(tours)
    tm = np.zeros((length, length))
    for i, tour in enumerate(tours):
        neighbours = list(get_neighbours(tour, move))
        n_neighbours = len(neighbours)
        for neighb in neighbours:
            j = tours.index(neighb)
            tm[i][j] += 1.0 / n_neighbours
    return tm

def get_tm_first_row(n, move="two_opt"):
    tours = list(tsp_tours(n))
    length = len(tours)
    tm = np.zeros((length,))
    tour = tours[0]
    neighbours = list(get_neighbours(tour, move))
    n_neighbours = len(neighbours)
    for neighb in neighbours:
        i = tours.index(neighb)
        tm[i] += 1.0 / n_neighbours
    return tm


def get_neighbours(t, move):
    """Iterate through all possible neighbours using the given type of move."""
    n = len(t)
    if move == "two_opt":
        for a in range(n):
            # b=a-1, b=a, and b=a+1 are invalid
            for b in range(a+2, a+n-1):
                b = b % n
                _a, _b = min(a, b), max(a, b)
                yield two_opt(t, (_a, _b))
    elif move == "twoh_opt":
        for a in range(n):
            b = (a+1) % n
            # c can be any node other than a or b
            for c in range(a+2, a+n):
                c = c % n
                yield twoh_opt(t, (a, b, c))
    elif move == "swap_two":
        for a in range(n):
            for b in range(a+1, n):
                yield swap_two(t, (a, b))
    elif move == "swap_adj":
        for a in range(n):
            b = (a+1) % n
            yield swap_adj(t, (a, b))
    elif move == "three_opt":
        for triple in _three_opt_choose_edges_iter(n):
            for which in [3, 4, 5, 6]:
                yield three_opt_deterministic(t, triple, which)
    elif move == "three_opt_broad":
        for triple in _three_opt_choose_edges_iter(n):
            for which in [1, 2, 3, 4, 5, 6, 7]:
                yield three_opt_deterministic(t, triple, which)
    else:
        raise ValueError("Unknown move " + move)


def sample_transitions(n, move="two_opt", nsamples=10000):
    length = len(list(tsp_tours(n)))
    tm = np.zeros((length, length))
    delta = 1.0 / nsamples

    if move == "three_opt":
        move = three_opt
    elif move == "three_opt_broad":
        move = three_opt_broad
    elif move == "two_opt":
        move = two_opt
    elif move == "swap_two":
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

def tsp_tm_wrapper(dirname, move="two_opt"):
    # dirname should be <dir>/tsp_length_6_2opt, for example
    t = dirname.find("tsp_length_")
    length = dirname[t:].split("_")[2]
    length = int(length)
    tm = get_tm(length, move)
    print(tm)
    outfilename = dirname + "/TP.dat"
    np.savetxt(outfilename, tm)
    kt = kendall_tau_permutation_distances(length)
    outfilename = dirname + "/KendallTau.dat"
    np.savetxt(outfilename, kt)

def count_permutations(n):
    # we canonicalise on starting at 0, so there are (n-1)! tours
    return math.factorial(n-1)

def test_op(op):
    n = 8
    m = 100000
    for i in range(m):

        p = range(n)
        p = op(p)
        yield tuple(p)



if __name__ == "__main__":
    # c = collections.Counter(test_op(three_opt))
    # print c.most_common(100)
    # print len(c.most_common(100))

    print(len(list(tsp_tours(6))))
    print(len(list(tsp_tours(7))))
    print(len(list(tsp_tours(8))))
