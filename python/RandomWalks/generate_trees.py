#!/usr/bin/env python

import sys
from itertools import product

vars = list("xy")
fns = {"+": 2, "-": 2, "*": 2, "/": 2}

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

if __name__ == "__main__":
    for item, d in trees_of_depth_LE(int(sys.argv[1])):
        print(item)
