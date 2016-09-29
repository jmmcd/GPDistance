import itertools
import random
from tsp import tsp_tours

"""A TSP problem is defined by a cost matrix C of size (n, n). A tour
is a permutation of the integers (0, n-1). The cost of a tour is the
sum of the appropriate n elements from C.

For a tour such as (0, 1, 2, 3, 4), we have an equation

C[0,1] + C[1, 2] + C[2, 3] + C[3, 4] + C[4, 0] = f(0, 1, 2, 3, 4)

Let's think of all elements of C reshaped as a row vector:

C[0,0] C[0,1] C[0,2] C[0,3] ... C[1,0] C[1,1] ... C[4,4]

Then we can write the LHS of the above equation as coefficients:

0 1 0 0 0   0 0 1 0 0  ...  1 0 0 0 0      |    f(0, 1, 2, 3, 4)

There are (n-1)! distinct tours, so we (n-1)! of these equations, and
there are n**2 values in C. So we get a binary matrix of shape
((n-1)!, n**2). Call this matrix a.

Suppose we are given a cost matrix C, we can calculate the costs of
all (n-1)!  tours. Call this vector b. It fully specifies a cost
function (assuming a fixed ordering on the tours, ie (01234),
(01243),...).

We can then solve the system aC = b for C, where we are again thinking
of C as a row vector. In this way, given only a vector b, we can
attempt to recover a C that would give rise to it.


Now, if TSP is *closed under permutation* (CUP) in the terms of the
sharpened No Free Lunch theorem, then we can permute the costs of
tours in any way we like, for example we can swap the costs of the
min-cost and max-cost tours, and the resulting matrix b1 will still
correspond to a TSP (potentially with a different cost matrix
C1). Thus we should be able to solve aC1 = b1 to recover a plausible
cost matrix. That would show that the vector b1 is a TSP problem.


This code is intended to explore this. We generate a random TSP cost
matrix, calculate the resulting b, construct the coefficient matrix a,
and then use Numpy.linalg.lstsq to solve the system. We then permute
to get b1, and again attempt to solve the system. The latter fails,
showing that the permuted b1 is not a TSP, hence TSP is not CUP.

However, a problem arises, which is that when solving the original
system aC = b, the solution that is found has some C < 0, which is not
the original C we created! To shore up this argument, I would like to
use an alternative method of solution, such as an LP solver or Sympy,
to impose the constraint that C > 0.

"""

# random TSP cost matrix
def random_tsp(n):
    return np.random.random((n, n))

# cost of a given tour, given a cost matrix
def cost(tour, C):
    i, j = tour_to_idx(tour)
    return C[i, j].sum()

# helper function: get the indices needed to pick out the cost matrix
# elements relevant to a given tour t
def tour_to_idx(t):
    ij = zip(t, t[1:] + [t[0]])
    return zip(*ij)

# helper function: get the indices needed to set the elements of the
# 'a' matrix
def tour_to_1d_idx(tour, n):
    return (i*n+j for i, j in zip(*tour_to_idx(tour)))

C = random_tsp(n)

# calculate cost of each tour
b = np.array([cost(t, x) for t in tsp_tours(n)])

# calculate permuted cost values
b1 = b.copy()
amx, amn = b1.argmax(), b1.argmin()
b1[amx], b1[amn] = b1.min(), b1.max()

# have a look at
for item in tour_to_1d_idx(t, n): print(item)

# construct a
a = np.zeros((b.shape[0], n*n))

for ti, tour in enumerate(tsp_tours(n)):
    a[ti, list(tour_to_1d_idx(tour, n))] = 1

# solve aC = b
C, res, rank, s = np.linalg.lstsq(a, b)
np.allclose(np.dot(a, C), b) # True, but notice C is not all positive!

# solve aC1 = b1
C1, res, rank, s = np.linalg.lstsq(a, b1)
np.allclose(np.dot(a, C1), b1) # False - no solution found




"""

Next idea is to try using sympy. We can re-use the constructions above,
and create X as a sympy MatrixSymbol.

"""

from sympy import Matrix
from sympy import MatrixSymbol


A = Matrix(a)
B = Matrix(b)

A.gauss_jordan_solve(B) # says "no solution", but it should have a solution

X = MatrixSymbol("X", 25, 1) # not sure about other approaches, like this one -- goal would be to impose the positive=True assumption on X
