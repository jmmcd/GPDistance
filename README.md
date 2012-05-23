GPDistance
==========

For measuring distances between genetic programming trees, both syntactically and via operators.

This is GPDistance, a set of tools in Java and Python for research
into distances between trees in genetic programming. It calculates two
types of distance: 

* Syntactic distances work by comparing the contents and structures of
  two trees. For now the syntactic distances are *tree-edit distance*
  (TED), *tree-alignment distance* (TAD), *overlap distance* (OVD),
  *feature-vector distance* (FVD), and *normalised compression
  distance* (NCD). GPDistance uses a Java implementation of the
  Zhang-Shasha tree-edit distance algorithm which is Copyright (C)
  2004 Stephen Wan (see
  http://web.science.mq.edu.au/~swan/howtos/treedistance/package.html).

* Operator distances work by studying the transition probabilities
  between pairs of nodes via particular operators. In this way we can
  calculate the *expected length of a random walk* (ELRW) and the
  *highest-probability path* (HPP). The ELRW is calculated with the
  help of ergodic.py, Copyright (C) 2012 Sergio J. Rey of PySAL (see
  http://code.google.com/p/pysal/), derived from the treatment in
  Kemeny, John, G. and J. Laurie Snell (1976) Finite Markov
  Chains. Springer-Verlag, Berlin.


Build
-----

There's a Makefile. Try $ make to compile everything. It calls into
Makefiles in sub-directories to compile them. The Python code requires
NumPy.

Usage
-----

See the Makefile to see some possibilities for what to run. Eg to run
a Metropolis-Hastings algorithm, try $ make mh. There are some scripts
for running multiple runs, as well as an R script for doing some
statistics on results.

TODO
----

There are a few assumptions in the code which would need to be checked
before extending it to consider (eg) other mutation operators, other
languages, etc.

Need a way to check how well-connected the graphs really are. When
using MH, it makes a *complete* graph over all generated nodes, so
everything is connected to everything else. But it would be
interesting to observe how many individuals are connected by the
actual mutations in the sampling -- ie how many of the
Metropolis-Hastings paths intersect -- before every pair of
individuals is connected up.


Publications
------------

If you wish to use and cite this work, please cite this earlier paper
which used many of the same concepts and methods (a newer publication
is in preparation):

McDermott, O'Reilly, Vanneschi, and Veeramachaneni, "How Far Is It
From Here to There? A Distance that is Coherent with GP Operators",
EuroGP 2011, Springer.

