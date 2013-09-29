#!/usr/bin/env python

import sys
from itertools import product

vars = list("xy")
fns = {"+": 2, "-": 2, "*": 2, "/": 2}

# Given this language, there are 2 trees of depth 0, 18 of depth 1,
# 1298 of depth 2, and 6739218 of depth 3. 

def trees_of_depth(n):
    """Generate all trees of exactly depth n, along with their
    depths."""
    if n == 0:
        for item in vars:
            yield item, 0
    else:
        for fn in fns:
            for children in product(trees_of_depth_LE(n-1), repeat=fns[fn]):
                if any(child[1] == n-1 for child in children):
                    yield ("(" + fn + " "
                           + " ".join(child[0] for child in children)
                           + ")"), n
                    
def trees_of_depth_LE(n):
    """Generate all trees up to and including depth n, along with
    their depths."""
    for d in range(n+1):
        for item, d1 in trees_of_depth(d):
            yield item, d1

def count_trees_of_depth_LE(n):
    """Count the number of trees of depth less than or equal to n."""
    return sum(count_trees_of_depth(i) for i in range(n+1))
            
def count_trees_of_depth(n):
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
        return len(fns.keys()) * (count_trees_of_depth(n-1) ** 2 + sum(
            2 * count_trees_of_depth(n-1) * count_trees_of_depth(m)
            for m in range(n-1)))

if __name__ == "__main__":
    n = int(sys.argv[1])
    if len(sys.argv) > 2 and sys.argv[2] == "enumerate":
        for item, d in trees_of_depth_LE(n):
            print(item)
    else:
        print(count_trees_of_depth_LE(n))
