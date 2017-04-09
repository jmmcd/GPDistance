"""Let's choose n = 3 variables, x0, x1, x2. So we also have y0, y1, y2 where yi is defined as NOT xi.

The search space X contains 2**3 = 8 points.

There are then 6 "variables", so (6 choose 2) = 15 possible clauses.

Any formula phi must consist of some subset of these clauses, so the number of possible formulae is 2**15 = 32768.

Another way of counting is: for m = 0 clauses, there is 1 possible fitness value (0); for m > 1 there are m+1 possible fitness values (0, 1, ... m); there may be some special cases eg for m = 15, phi will include all possible clauses, and some will be incompatible, so I guess f(x) = 15 is never achieved. Enumerate these to check. But for a first pass, say that for m clauses, we have m+1 possible fitness values, hence 8**(m+1) possible fitness functions.

For m clauses, we have (15 choose m) possible formulae.

Let's save all the formulae and all the tables that we achieve

And also create the dictionaries mapping formula <-> fitness table, eg

[(x0, y1), (x1, y2)] <-> [0, 1, 0, 2, 2, 2, 0, 1]

Then for each table, we try all its permutations until we find one that is not reachable using any formula.
"""

from itertools import product, combinations, permutations
import scipy.misc
choose = lambda N, k: int(scipy.misc.comb(N, k))

n = 3 # number of variables
varnames = ["x[%d]" % i for i in range(n)] + ["y[%d]" % i for i in range(n)]

def phi

def make_fit(phi):
    """phi will be a string of the form "[(x[0], y[1]), (x[1], y[2])]"

    >>> make_fit("[(x[0], y[1]), (x[1], y[2])]")([0, 0, 0])
    2
    >>> make_fit("[(x[0], y[1]), (x[1], y[2])]")([1, 1, 1])
    2
    """

    def fitness(x):
        y = [not xi for xi in x]
        phi_ = eval(phi)
        return sum((a or b) for a, b in phi_)

    return fitness

f = make_fit("[(x[0], y[1]), (x[1], y[2])]")

def fit_table(phi):
    f = make_fit(phi)
    return tuple([f(x) for x in product([0, 1], repeat=n)])

def enum_formulae(n, m=None):
    if m is None:
        nclauses = choose(2 * n, 2)
        for m in range(1, nclauses + 1):
            print(m)
            yield from enum_formulae(n, m)
    else:
        clauses = ["(%s, %s)" % (a, b) for a, b in combinations(varnames, 2)]
        formulae = ["[" + ", ".join(comb) + "]" for comb in combinations(clauses, m)]
        yield from formulae


def search_for_unreachable(fit_tables, fit_tables_to_formulae):
    for i, table in enumerate(fit_tables):
        print(str(i) + " of " + str(len(fit_tables)))
        print(table)
        print(fit_tables_to_formulae[table])
        for perm in permutations(table):
            if perm not in fit_tables:
                yield table, perm


def main():
    fit_tables = []
    formulae = []
    fit_tables_to_formulae = dict()
    formulae_to_fit_tables = dict()
    nformulae = 0
    for phi in enum_formulae(n):
        #f = make_fit(phi)
        table = fit_table(phi)
        #print(phi + ": " + str(table))
        fit_tables.append(table)
        formulae.append(phi)
        fit_tables_to_formulae[table] = phi
        formulae_to_fit_tables[phi] = table
        nformulae += 1

    print(len(fit_tables_to_formulae)) # 9922
    print(nformulae) # 32767

    for table, perm in search_for_unreachable(fit_tables, fit_tables_to_formulae):
        print("Table " + str(table) + " is reachable")
        print("by formula " + str(fit_tables_to_formulae[table]))
        print("but its permutation " + str(perm) + " is not")
        raise StopIteration
