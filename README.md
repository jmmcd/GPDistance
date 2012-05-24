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
  distance* (NCD). The Java implementation of the Zhang-Shasha
  tree-edit distance algorithm is Copyright (C) 2004 Stephen Wan (see
  http://web.science.mq.edu.au/~swan/howtos/treedistance/package.html). It
  also uses a Java transliteration of the Clojure overlap distance by
  Lee Spector (see http://push.i3ci.hampshire.edu/2011/06/07/overlap/).

* Operator distances work by studying the transition probabilities
  between pairs of nodes via particular operators. In this way we can
  calculate the *expected length of a random walk* (ELRW), the
  *highest-probability path* (HPP), and the *shortest path length*
  (SPL). The ELRW is calculated with the help of ergodic.py, Copyright
  (C) 2012 Sergio J. Rey of PySAL (see
  http://code.google.com/p/pysal/), derived from the treatment in
  Kemeny, John, G. and J. Laurie Snell (1976) Finite Markov
  Chains. Springer-Verlag, Berlin.


Build
-----

There's a Makefile. Try `make` to compile everything. It calls into
Makefiles in sub-directories to compile them. The Python code requires
NumPy and SciPy. Use `apt` or your favourite package manager, or try
`pip install numpy scipy`, but be warned that installing them is known
to cause pain especially on OSX.


Usage
-----

See the root Makefile to see some possibilities for what to run. Eg to
write matrices with all pairwise distances for the complete space of
depth 1, try `make completeMatrices1`.


TODO
----

There are a few assumptions in the code which would need to be checked
before extending it to consider (eg) other mutation operators, other
languages, etc.


Publications
------------

If you wish to use and cite this work, please cite this earlier paper
which used many of the same concepts and methods (a newer publication
is in preparation):

McDermott, O'Reilly, Vanneschi, and Veeramachaneni, "How Far Is It
From Here to There? A Distance that is Coherent with GP Operators",
*EuroGP 2011*, Springer.

