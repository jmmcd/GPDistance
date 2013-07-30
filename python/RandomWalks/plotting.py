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
    gold_names = ["D_TP", "MFPT", "SP", "STEPS"]
    
    for name in gold_names + syntactic_distance_names:
        w = np.genfromtxt(codename + "/" + name + ".dat")
        assert(len(w) == len(names))
        make_grid(w, False, codename + "/" + name)
    print names # better to print them in a list somewhere than in the graph
        
def make_grid(w, names, filename):
    # TODO for now we dont rescale the data, although there are some
    # matrices which would benefit from a log transform or
    # similar. matshow() internally scales the data so that the
    # smallest numbers go to black and largest to white.

    fig = plt.figure(figsize=(25, 25))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_frame_on(False)
    im = ax.matshow(w, cmap=cm.gray, interpolation="none")

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
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png")

def get_kendall_tau(x, y):
    """Return Kendall's tau, a non-parametric test of association. If
     one of the variables is constant, a FloatingPointError will
     happen and we can just say that there was no association. Note
     this runs Kendall's tau-b, accounting for ties and suitable for
     square tables:
     [http://en.wikipedia.org/wiki/Kendall_tau_rank_correlation_coefficient#Tau-b]
     [http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kendalltau.html]"""

    if "scipy" not in locals():
        import scipy.stats
    # make sure we raise any error (so we can catch it), don't just
    # splat it on the terminal
    old = np.seterr(all='raise')
    try:
        corr, p = scipy.stats.kendalltau(x, y)
    except (FloatingPointError, TypeError):
        # TypeError can arise in Kendall tau when both x and y are
        # constant
        corr, p = 0.0, 1.0
    # restore old error settings
    np.seterr(**old)
    return corr, p

def get_spearman_rho(x, y):
    """Return Spearman's rho, a non-parametric test of association. If
     one of the variables is constant, a FloatingPointError will
     happen and we can just say that there was no association.
     [http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient].
     The reason for using this is that the usual Pearson's correlation
     assumes normal distributions, which our distances certainly
     aren't, and Kendall's tau is O(n^2), so unbearably slow."""

    if "scipy" not in locals():
        import scipy.stats
    # Make sure we raise any error (so we can catch it), don't just
    # splat it on the terminal. However we will ignore underflow
    # because it seems to happen with every single one. Similar
    # approach here:
    # [https://code.google.com/p/wnd-charm/source/browse/pychrm/trunk/pychrm/FeatureSet.py?r=723]
    old = np.seterr(all='raise', under='ignore')
    try:
        corr, p = scipy.stats.spearmanr(x, y)
    except FloatingPointError:
        corr, p = 0.0, 1.0
    # restore old error settings
    np.seterr(**old)
    return corr, p

def get_pearson_r(x, y):
    """Return Pearson's R, the most common test of correlation, which
     assumes the variables are each normally distributed. If one of
     the variables is constant, a FloatingPointError will happen and
     we can just say that there was no association."""

    if "scipy" not in locals():
        import scipy.stats
    # Make sure we raise any error (so we can catch it), don't just
    # splat it on the terminal. 
    old = np.seterr(all='raise')
    try:
        corr, p = scipy.stats.pearsonr(x, y)
    except FloatingPointError:
        corr, p = 0.0, 1.0
    # restore old error settings
    np.seterr(**old)
    return corr, p

def make_correlation_tables(codename, txt=""):

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

    gold_names = ["D_TP", "MFPT", "SP", "STEPS"]
    d = {}
    for name in syntactic_distance_names + gold_names:
        # print("reading " + name)
        m = np.genfromtxt(codename + "/" + name + ".dat")
        d[name] = m.reshape(len(m)**2)

    def do_line(syn):
        line = syn.replace("_TP", r"$_{\mathrm{TP}}$")
        for gold in gold_names:
            print("getting association between " + gold + " " + syn)
            if corr_type == "spearmanrho":
                corr, p = get_spearman_rho(d[gold], d[syn])
            elif corr_type == "kendalltau":
                corr, p = get_kendall_tau(d[gold], d[syn])
            elif corr_type == "pearsonr":
                corr, p = get_pearson_r(d[gold], d[syn])
            else:
                print("Unknown correlation type " + corr_type)
            if p < 0.05:
                sig = "*"
            else:
                sig = " "
            line += r" & {0:1.2f} \hfill {1} ".format(corr, sig)
        line += r"\\"
        f.write(line + "\n")
        
    for corr_type in ["spearmanrho", "pearsonr", "kendalltau"]:
        if corr_type == "kendalltau" and len(d["STEPS"]) > 1000:
            print("Omitting Kendall tau because it is infeasible for large matrices")
            continue
        filename = codename + "/correlation_table_" + corr_type + ".tex"
        f = open(filename, "w")

        f.write(r"""\begin{table}
\centering
\caption{Correlations between distance measures """)
        if corr_type == "spearmanrho":
            f.write("using Spearman's rho: ")
        elif corr_type == "pearsonr":
            f.write("using Pearson's R: ")
        elif corr_type == "kendalltau":
            f.write("using Kendall's tau: ")
        f.write(txt)
        f.write(r"""\label{tab:correlationresults_""" + os.path.basename(codename) + "}}\n")
        f.write(r"""\begin{tabular}{l|l|l|l|l}
     & D$_{\mathrm{TP}}$ & MFPT & SP & STEPS \\
\hline
\hline
""")

        for syn in syntactic_distance_names:
            do_line(syn)
        f.write(r"""\hline
     \hline
""")
        for gold in gold_names:
            do_line(gold)

        f.write(r"""\end{tabular}
    \end{table}""")

    f.close()

def compare_sampled_v_calculated(codename):
    if "scipy" not in locals():
        import scipy.stats
    
    stp = np.genfromtxt(codename + "/TP_sampled.dat")
    stp = stp.reshape(len(stp)**2)
    tp = np.genfromtxt(codename + "/TP.dat")
    tp = tp.reshape(len(tp)**2)

    filename = codename + "/compare_TP_calculated_v_sampled.tex"
    f = open(filename, "w")
    corr, p = scipy.stats.pearsonr(tp, stp)
    f.write("Pearson R correlation " + str(corr) + "; ")
    f.write("p-value " + str(p) + ". ")
    corr, p = scipy.stats.spearmanr(tp, stp)
    f.write("Spearman rho correlation " + str(corr) + "; ")
    f.write("p-value " + str(p) + ". ")
    if len(stp) < 1000:
        corr, p = scipy.stats.kendalltau(tp, stp)
        f.write("Kendall tau correlation " + str(corr) + "; ")
        f.write("p-value " + str(p) + ". ")
    else:
        f.write("Omitting Kendall tau because it is infeasible for large matrices. ")
    f.close()

def write_steady_state(codename):
    """Read in a TP matrix given a codename. Use ergodic.steady_state
    to calculate the long-run steady-state, which is a vector
    representing how long the system will spend in each state in the
    long run. If not uniform, that is a bias imposed by the operator
    on the system. Write it out and plot it. Calculate the stddev of
    the steady-state as well."""
    import ergodic
    names = open(codename + "/all_trees.dat").read().strip().split("\n")
    names = map(lambda x: x.strip(), names)
    tp = np.genfromtxt(codename + "/TP.dat")
    ss = ergodic.steady_state(np.matrix(tp))
    ss = np.real(ss) # discard zero imaginary parts
    print("Stddev of steady-state for " + codename + ": ")
    print(np.std(ss))
    print("Sum of steady-state for " + codename + ": ")
    print(np.sum(ss))
    fig = plt.figure(figsize=(6.5, 3.5))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    offset = int(log(len(ss), 2)) + 1
    ax.set_xlim((-offset, len(ss)-1+offset))
    plt.ylabel("Log(steady-state probability)")
    ax.plot(ss)
    ax.set_xticklabels([], [])
    filename = codename + "/steady_state.dat"
    np.savetxt(filename, ss)
    filename = codename + "/steady_state"
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png")
    

if __name__ == "__main__":
    cmd = sys.argv[1]
    codename = sys.argv[2]

    if cmd == "compareTPCalculatedVSampled":
        compare_sampled_v_calculated(codename)
    elif cmd == "makeCorrelationTable":
        txt = sys.argv[3]
        make_correlation_tables(codename, txt)
    elif cmd == "makeGridPlots": 
        make_grid_plots(codename)
    elif cmd == "writeSteadyState":
        write_steady_state(codename)
