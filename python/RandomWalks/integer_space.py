#!/usr/bin/env python

import numpy as np
import random
import itertools

def int_space_make_ops(N):
    def uniform(x):
        return random.randrange(N)
    def ten_percent(x):
        width = N / 10
        choice = random.randrange(-width, width)
        return (x + choice) % N
    def deterministic(x):
        return (x+1) % N
    return uniform, ten_percent, deterministic

def int_space_make_rows(N):
    x_unif = np.ones(N) / float(N)

    x_10pc = np.zeros(N)
    width = N / 10
    x_10pc[0:width] = 1.0 / width

    x_det = np.zeros(N)
    x_det[0] = 1.0

    return x_unif, x_10pc, x_det



# class IntegerSpace:
#     def __init__(self, n, op):
#         self.n = n
#         self.N = n
#         self.op = op

#     def mutate(self, x):
#         return self.op(x)

#     def write_TP_row0(self):
#         pass

#     def fitness(self, x):
#         return
