#!/usr/bin/env python2.7
# 2.7 because scipy doesn't import on other versions, on my imac

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.ticker import FuncFormatter, MaxNLocator, IndexLocator
import sys
import os
from math import *
import random

# MAXTICKS is 1000 in IndexLocator
class MyLocator(mpl.ticker.IndexLocator):
    MAXTICKS=1500

def make_grid_plots(codename):
    names = open(codename + "/all_trees.dat").read().strip().split("\n")
    names = map(lambda x: x.strip(), names)
    
    # gp distances
    syntactic_distance_names = [
        "NCD", "FVD",
        "NodeCount", "MinDepth", "MeanDepth", "MaxDepth",
        "Symmetry", "MeanFanout", "DiscreteMetric",
        "TED",
        "TAD0", "TAD1", "TAD2", "TAD3", "TAD4", "TAD5",
        "OVD"
    ]
    # ga distances
    if "ga_length" in codename:
        syntactic_distance_names = [
            "Hamming"
        ]
    gold_names = ["D_TP", "FMPT", "SP", "STEPS"]
    
    for name in gold_names + syntactic_distance_names:
        w = np.genfromtxt(codename + "/" + name + ".dat")
        assert(len(w) == len(names))
        make_grid(w, names, codename + "/" + name + ".pdf")

def make_grid(w, names, filename):
    # rescale/reshape: avoid uniform colours, and make higher numbers
    # darker
    de = -np.log(w)

    fig = plt.figure(figsize=(25, 25))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_frame_on(False)
    im = ax.matshow(de, cmap=cm.gray, interpolation="none")

    if names:
        # Turn labels on
        ax.xaxis.set_major_locator(MyLocator(1, 0))
        ax.yaxis.set_major_locator(MyLocator(1, 0))
        ax.set_xticklabels(names, range(len(names)), rotation=90, size=1.0)
        ax.set_yticklabels(names, range(len(names)), size=1.0)
    else:
        # Turn them off
        ax.set_xticklabels([], [])
        ax.set_yticklabels([], [])
    
    ax.tick_params(length=0, pad=3.0)
    fig.savefig(filename)

def get_kendall_tau(x, y):
    """Return Kendall's tau, a non-parametric test of association. If
     one of the variables is constant, a FloatingPointError will
     happen and we can just say that there was no association. Note
     this runs Kendall's tau-b, accounting for ties and suitable for
     square tables:
     [http://en.wikipedia.org/wiki/Kendall_tau_rank_correlation_coefficient#Tau-b]
     [http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kendalltau.html]

    """

    if "scipy" not in locals():
        import scipy.stats
    # make sure we raise any error (so we can catch it), don't just
    # splat it on the terminal
    old = np.seterr(all='raise')
    try:
        corr, p = scipy.stats.mstats.kendalltau(x, y)
    except FloatingPointError:
        corr, p = 0.0, 1.0
    # restore old error settings
    np.seterr(**old)
    return corr, p

def make_correlation_table(codename, txt=""):
    # gp distances
    syntactic_distance_names = [
        "NCD", "FVD",
        "NodeCount", "MinDepth", "MeanDepth", "MaxDepth",
        "Symmetry", "MeanFanout", "DiscreteMetric",
        "TED",
        "TAD0", "TAD1", "TAD2", "TAD3", "TAD4", "TAD5",
        "OVD"
    ]
    # ga distances
    if "ga_length" in codename:
        syntactic_distance_names = [
            "Hamming"
        ]

    print(r"""\begin{table}
\centering
\caption{Correlations between distance measures: """ + txt)
    print(r"""\label{tab:correlationresults_""" + os.path.basename(codename) + r"}}")
    print(r"""\begin{tabular}{l|l|l|l|l}
 & D$_{\mathrm{TP}}$ & FMPT & SP & STEPS \\
\hline
\hline""")

    gold_names = ["D_TP", "FMPT", "SP", "STEPS"]
    d = {}
    array_len = -1
    for name in syntactic_distance_names + gold_names:
        print("reading " + name)
        d[name] = np.genfromtxt(codename + "/" + name + ".dat")
        array_len = len(d[name])

    # sample_positions = random.sample([(i, j)
    #                                   for i in range(array_len)
    #                                   for j in range(array_len)], 1000)
    # def sample_from(d):
    #     return np.array([d[i][j] for i, j in sample_positions])

    def do_line(syn):
        line = syn.replace("_TP", r"$_{\mathrm{TP}}$")
        for gold in gold_names:
            corr, p = get_kendall_tau(d[gold], d[syn])
            if p < 0.05:
                sig = "*"
            else:
                sig = " "
            line += r" & {0:1.2f} \hfill {1} ".format(corr, sig)
        line += r"\\"
        print(line)
        
    for syn in syntactic_distance_names:
        do_line(syn)
    print(r"""\hline
 \hline""")
    for gold in gold_names:
        do_line(gold)

    print(r"""\end{tabular}
\end{table}""")

if __name__ == "__main__":
    codename = sys.argv[1]
    
    # txt = sys.argv[2]
    # make_correlation_table(codename, txt)
    
    make_grid_plots(codename)
