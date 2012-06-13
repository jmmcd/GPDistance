#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.ticker import FuncFormatter, MaxNLocator, IndexLocator
import scipy.stats
import sys
import os
from math import *
np.seterr(all='raise')

# MAXTICKS is 1000 in IndexLocator
class MyLocator(mpl.ticker.IndexLocator):
    MAXTICKS=1500

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
    for name in syntactic_distance_names + gold_names:
        d[name] = np.genfromtxt(codename + "/" + name + ".dat")

    def do_line(syn):
        line = syn.replace("_TP", r"$_{\mathrm{TP}}$")
        for gold in gold_names:
            try:
                corr, p = scipy.stats.mstats.kendalltau(d[gold], d[syn])
            except FloatingPointError:
                corr, p = 0.0, 1.0
            # corr = abs(corr)
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
    # for fast testing of aesthetic changes
    # w = np.array([[3, 4, 5, 6], [1, 2, 3, 4], [2, 3, 4, 5], [6, 7, 8, 9]])

    codename = sys.argv[1]
    txt = sys.argv[2]
    make_correlation_table(codename, txt)

    
    # names = open(codename + "_trees.txt").read().strip().split("\n")
    # names = map(lambda x: x.strip(), names)

    # #w = np.genfromtxt(codename + "_tm.dat")
    # for suffix in ["mean", "len", "min", "max", "std"]:
    #     w = np.genfromtxt("randomwalking_x{0}.dat".format(suffix))
    #     assert(len(w) == len(names))
    #     #make_grid(w, None, codename + "_tm.pdf")
    #     make_grid(w, None, "randomwalking_x{0}.pdf".format(suffix))