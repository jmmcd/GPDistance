#!/usr/bin/env python

import sys
from itertools import product
from collections import OrderedDict
import numpy as np
from numpy import add, subtract, multiply, divide, sin, cos, exp, log, power, square, sqrt
np.seterr(all='raise')

def trees_of_depth(n, vars, fns, as_string=True):
    """Generate all trees of exactly depth n, along with their
    depths."""
    if n == 0:
        for item in vars:
            yield item, 0
    else:
        for fn in fns:
            for children in product(
                trees_of_depth_LE(n-1, vars, fns, as_string), repeat=fns[fn]):
                if any(child[1] == n-1 for child in children):
                    if as_string:
                        yield ("(" + fn + " "
                               + " ".join(child[0] for child in children)
                               + ")"), n
                    else:
                        yield [fn] + [child[0] for child in children], n

def trees_of_depth_LE(n, vars, fns, as_string=True):
    """Generate all trees up to and including depth n, along with
    their depths."""
    for d in range(n+1):
        for item, d1 in trees_of_depth(d, vars, fns, as_string):
            yield item, d1

def count_trees_of_depth_LE(n, vars, fns):
    """Count the number of trees of depth less than or equal to n."""
    return sum(count_trees_of_depth(i, vars, fns)
               for i in range(n+1))

def count_trees_of_depth(n, vars, fns):
    """Count the trees of depth exactly n."""
    if n == 0:
        return len(vars)
    else:
        # there are len(fns.keys()) ways to choose the root

        # then you can choose trees of depth n-1 for *both* subtrees.
        # that gives count_trees_of_depth(n-1) ** 2

        # otherwise, you must choose a left subtree of depth n-1, and
        # a right subtree of any depth less than that. that gives the
        # sum of count_trees_of_depth(n-1) * count_trees_of_depth(m),
        # with m ranging from 0 up to n-2 inclusive.

        # or you can do vice versa (right subtree of depth n-1, left
        # subtree of any depth less than that). so multiply previous
        # quantity by 2.
        return len(fns.keys()) * (
            count_trees_of_depth(n-1, vars, fns) **
            2 + sum(2 * count_trees_of_depth(n-1, vars, fns)
                    * count_trees_of_depth(m, vars, fns)
                    for m in range(n-1)))

def count_trees_of_given_shape(t, vars, fns):
    """Tree t is a string, eg (* (+ x y) y). This has n = 2 internal
    nodes, m = 3 leaves. Number of trees of the same shape will be
    4^n2^m (if there are 4 possibilities for internal nodes and 2 for
    leaves)."""
    def noccurrences(s1, s2):
        """How many times do items in s2 occur in s1?"""
        return sum(s1.count(a) for a in s2)
    internal = vars
    external = "".join(fns.keys())
    n = noccurrences(t, internal)
    m = noccurrences(t, external)
    return (len(internal)**n) * (len(external)**m)

def shapes_of_depth_LE(n):
    for t in trees_of_depth_LE(n, "x", {"*": 2}): yield t

# AQ is the analytic quotient from: Ji Ni and Russ H. Drieberg and
# Peter I. Rockett, "The Use of an Analytic Quotient Operator in
# Genetic Programming", IEEE Transactions on Evolutionary Computation
def AQ(x, y):
    return x/sqrt(1.0+y*y)

def semantics(t):
    f = tree_to_fn(t)
    x = np.linspace(-5.0, 5.0, 26)
    X = np.meshgrid(x, x)
    X = np.array([X[0].ravel(), X[1].ravel()])
    return X, f(X)
    
def fitness(t):
    X, t_vals = semantics(t)
    target = target_fn(X)
    return np.sqrt(np.mean((target - t_vals)**2))
    
def target_fn(X):
    """Pagie and Hogeweg"""
    return 1.0 / (1.0 + X[0]**-4.0) + 1.0 / (1.0 + X[1]**-4.0)
    
def tree_to_fn(t):
    if type(t) == str:
        # it's either X[0] or similar or a const: just eval
        if t == "x" or t == "x0":
            t = "X[0]"
        elif t == "y" or t == "x1":
            t = "X[1]"
        f = lambda X: eval(t)
    else:
        # it's a fn -- make fn, eval args, then apply
        op = str(t[0])
        if op == "+": op = add
        elif op == "*": op = multiply
        elif op == "/": op = AQ
        elif op == "-": op = subtract
        f = lambda X: op(*(tree_to_fn(ti)(X) for ti in t[1:]))
    return f

def enumerate_semantics(n, vars, fns):
    ntrees = count_trees_of_depth_LE(2, vars, fns)
    result = []
    for t, d in trees_of_depth_LE(n, vars, fns, False):
        result.append(semantics(t)[1])
    result = np.array(result)
    return result

def semantic_distances(n, vars, fns):
    semantics = enumerate_semantics(n, vars, fns)
    ntrees = semantics.shape[0]
    print "ntrees", ntrees
    result = np.zeros((ntrees, ntrees))
    for i, t in enumerate(semantics):
        print "len t", len(t)
        for j, s in enumerate(semantics):
            d = np.sqrt(np.mean((t - s)**2.0))
            result[i, j] = d
    return result
    
if __name__ == "__main__":

    # Given this language, there are 2 trees of depth 0, 18 of depth
    # 1, 1298 of depth 2, and 6739218 of depth 3.
    vars = ["x0", "x1"]
    fns = OrderedDict([("*", 2), ("+", 2), ("-", 2), ("/", 2)])

    n = int(sys.argv[1])
    if len(sys.argv) > 2 and sys.argv[2] == "enumerate":
        for item, d in trees_of_depth_LE(n, vars, fns, False):
            print(item)
            
    elif len(sys.argv) > 2 and sys.argv[2] == "shapes":
        for item, d in shapes_of_depth_LE(n):
            print("%d: %s" % (count_trees_of_given_shape(item, vars, fns),
                             item))
            
    elif len(sys.argv) > 2 and sys.argv[2] == "enumerate_semantics":
        dirname = sys.argv[3]
        result = enumerate_semantics(n, vars, fns)
        np.savetxt(dirname + "/all_trees_semantics.dat", result)
        
    elif len(sys.argv) > 2 and sys.argv[2] == "semantic_distances":
        dirname = sys.argv[3]
        result = semantic_distances(n, vars, fns)
        np.savetxt(dirname + "/SEMD.dat", result)
        
    else:
        print(count_trees_of_depth_LE(n, vars, fns))
