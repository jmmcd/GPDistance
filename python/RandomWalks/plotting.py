#!/usr/bin/env python2.7
# 2.7 because scipy doesn't import on other versions, on my imac

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.manifold import MDS
import numpy as np
from matplotlib.ticker import FuncFormatter, MaxNLocator, IndexLocator
import sys
import os
from math import *
import random

# MAXTICKS is 1000 in IndexLocator
class MyLocator(mpl.ticker.IndexLocator):
    MAXTICKS=1500

graph_distance_names = ["D_TP", "SD_TP",
                        "MFPT", "CT", "MFPT_VLA", "CT_VLA",
                        "SP", "STEPS",
                        ]
graph_distance_tex_names = ["$D_{\mathrm{TP}}$", r"S$D_{\mathrm{TP}}$",
                            "MFPT", "CT", "MFPT-VLA", "CT-VLA",
                            "SP", "STEPS",
                            ]

def make_grid_plots(dirname):
    if "depth" in dirname:
        names = open(dirname + "/all_trees.dat").read().strip().split("\n")
        names = map(lambda x: x.strip(), names)
    else:
        # Assume GA
        length = int(dirname.strip("/").split("_")[2])
        names = [bin(i)[2:] for i in range(2**length)]
    
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
    if "ga_length" in dirname:
        syntactic_distance_names = [
            "Hamming"
        ]
    
    for name in graph_distance_names + syntactic_distance_names:
        # if name not in ("OVD", "TAD0"): continue
        w = np.genfromtxt(dirname + "/" + name + ".dat")
        assert(len(w) == len(names))
        make_grid(w, False, dirname + "/" + name)
    print names # better to print them in a list somewhere than in the graph
        
def make_grid(w, names, filename):
    # TODO for now we dont rescale the data, although there are some
    # matrices which would benefit from a log transform or
    # similar. matshow() internally scales the data so that the
    # smallest numbers go to black and largest to white.

    # A uniform array will cause a "RuntimeWarning: invalid value
    # encountered in divide" when calculating the colorbar. So ignore
    # that.
    old = np.seterr(invalid='ignore')

    # Can put NaN on the diagonal to avoid plotting it -- makes a bit
    # more space available for other data. But misleading.
    
    # w += np.diag(np.ones(len(w)) * np.nan)

    figsize = (10, 10)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    # ax.set_frame_on(False)
    # consider other colour maps: cm.gray_r for reversed, autumn, hot,
    # gist_earth, copper, ocean, some others, or a custom one for
    # nicer images (not for publication, maybe). 
    im = ax.matshow(w, cmap=cm.gray, interpolation="nearest")
    fig.colorbar(im, shrink=0.775)
    
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
    fig.savefig(filename + ".pdf", dpi=100)
    fig.savefig(filename + ".png", dpi=100)

    # restore old error settings
    np.seterr(**old)


def get_kendall_tau(x, y):
    """Return Kendall's tau, a non-parametric test of association. If
     one of the variables is constant, a FloatingPointError will
     happen and we can just say that there was no association. Note
     this runs Kendall's tau-b, accounting for ties and suitable for
     square tables:
     [http://en.wikipedia.org/wiki/Kendall_tau_rank_correlation_coefficient#Tau-b]
     [http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kendalltau.html].
     However it is of n^2 time complexity, hence unsuitable for even
     medium-sized inputs."""

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

def make_correlation_tables(dirname, txt=""):

    # gp distances
    syntactic_distance_names = [
        "NCD", "FVD",
        "NodeCount", "MinDepth", "MeanDepth", "MaxDepth",
        "Symmetry", "MeanFanout", "DiscreteMetric",
        "TED",
        "TAD0", "TAD1", "TAD2", "TAD3", "TAD4", "TAD5",
        "OVD",
    ]
    # ga distances
    if "ga_length" in dirname:
        syntactic_distance_names = [
            "Hamming"
        ]

    d = {}
    for name in syntactic_distance_names + graph_distance_names:
        # print("reading " + name)
        m = np.genfromtxt(dirname + "/" + name + ".dat")
        d[name] = m.reshape(len(m)**2)

    def do_line(dist, dist_name):
        line = dist_name.replace("_TP", r"$_{\mathrm{TP}}$")
        for graph_distance in graph_distance_names:
            print("getting association between " + graph_distance + " " + dist)
            if corr_type == "spearmanrho":
                corr, p = get_spearman_rho(d[graph_distance], d[dist])
            elif corr_type == "kendalltau":
                corr, p = get_kendall_tau(d[graph_distance], d[dist])
            elif corr_type == "pearsonr":
                corr, p = get_pearson_r(d[graph_distance], d[dist])
            else:
                print("Unknown correlation type " + corr_type)
            if corr > 0.0 and p < 0.05:
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
        filename = dirname + "/correlation_table_" + corr_type + ".tex"
        f = open(filename, "w")

        f.write(r"""\begin{table*}
\centering
\caption{Correlations between distance measures """)
        if corr_type == "spearmanrho":
            f.write("using Spearman's rho: ")
        elif corr_type == "pearsonr":
            f.write("using Pearson's R: ")
        elif corr_type == "kendalltau":
            f.write("using Kendall's tau: ")
        f.write(txt)
        f.write(r"""\label{tab:correlationresults_""" + os.path.basename(dirname) + "}}\n")
        f.write(r"""\begin{tabular}{""" + "|".join("l" for i in range(1+len(graph_distance_names))) + r"""}
     & """ + " & ".join(graph_distance_tex_names)  + r"""\\
\hline
\hline
""")

        for name, tex_name in zip(graph_distance_names, graph_distance_tex_names):
            do_line(name, tex_name)
        f.write(r"""\hline
\hline
""")
        for syn in syntactic_distance_names:
            do_line(syn, syn)

        f.write(r"""\end{tabular}
\end{table*}""")

    f.close()

def compare_sampled_v_calculated(dirname):
    if "scipy" not in locals():
        import scipy.stats
    
    stp = np.genfromtxt(dirname + "/TP_sampled.dat")
    stp = stp.reshape(len(stp)**2)
    tp = np.genfromtxt(dirname + "/TP.dat")
    tp = tp.reshape(len(tp)**2)

    filename = dirname + "/compare_TP_calculated_v_sampled.tex"
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

def write_steady_state(dirname):
    """Get steady state, write it out and plot it. Also calculate the
    in-degree of each node, by summing columns. Plot that on the same
    plot. Calculate the correlation between steady-state and
    in-degree. Calculate the stddev of the steady-state as well, and
    (why not) the stddev of the TP matrix as well."""
    from random_walks import get_steady_state
    if "scipy" not in locals():
        import scipy.stats
    tp = np.genfromtxt(dirname + "/TP.dat")
    ss = get_steady_state(tp)
    s = ("Stddev " + str(np.std(ss)) + ". ")
    open(dirname + "/steady_state.tex", "w").write(s)
    s = ("Stddev " + str(np.std(tp)) + ". ")
    open(dirname + "/tp_stddev.tex", "w").write(s)
    cs = np.sum(tp, axis=0)
    cs /= np.sum(cs)
    s = ("Pearson correlation between steady-state vector "
         + "and normalised in-degree vector"
         + str(scipy.stats.pearsonr(ss, cs)) + ". ")
    open(dirname + "/in_degree.tex", "w").write(s)
    
    fig = plt.figure(figsize=(6.5, 3.5))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    offset = int(log(len(ss), 2)) + 1
    ax.set_xlim((-offset, len(ss)-1+offset))
    plt.ylabel("Log-probability")
    ax.plot(ss, label="Steady-state", lw=2)
    ax.plot(cs, label="Normalised in-degree", lw=2)
    ax.set_xticklabels([], [])
    if ("depth_1" in dirname) or ("ga_length" in dirname):
        plt.legend(loc=3)
    else:
        plt.legend(loc=1)
    filename = dirname + "/steady_state.dat"
    np.savetxt(filename, ss)
    filename = dirname + "/in_degree.dat"
    np.savetxt(filename, cs)
    filename = dirname + "/steady_state"
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png")

def make_mds_images(dirname):
    """Make MDS images for multiple distance matrices. Each matrix
    must be symmetric. Must not contain any infinities, which prevents
    SD_TP in a space like ga_length_4_per_ind."""

    for name in ["CT", "SD_TP", "TED", "TAD0", "OVD", "FVD"]:
        m = np.genfromtxt(dirname + "/" + name + ".dat")
        make_mds_image(m, dirname + "/" + name + "_MDS")
    
def make_mds_image(m, filename):
    """Given a matrix of distances, project into 2D space using
    multi-dimensional scaling and produce an image."""

    # Construct MDS object with various defaults including 2d
    mds = MDS(dissimilarity="precomputed")
    # Fit
    try:
        f = mds.fit(m)
    except ValueError:
        print("Can't run MDS for " + filename + " because it contains infinities.")
        return
    
    # Get the embedding in 2d space
    p = f.embedding_

    # Make an image
    plt.figure()
    # x- and y-coordinates
    plt.axes().set_aspect('equal')
    plt.scatter(p[:,0], p[:,1],
                marker='.', c='b')
    # For triangles and random colours use this:
    # marker='^',
    # c=[random.random() for i in range(len(p[:,0]))],
    # cmap=cm.autumn)
    plt.savefig(filename + ".png")
    plt.savefig(filename + ".pdf")

    # Could use custom marker types and colours to get an interesting
    # diagram?
    # http://matplotlib.org/api/path_api.html#matplotlib.path.Path
    # (pass marker=Path() from above), and then use random colours?

if __name__ == "__main__":
    cmd = sys.argv[1]
    dirname = sys.argv[2]

    if cmd == "compareTPCalculatedVSampled":
        compare_sampled_v_calculated(dirname)
    elif cmd == "makeCorrelationTable":
        txt = sys.argv[3]
        make_correlation_tables(dirname, txt)
    elif cmd == "makeGridPlots": 
        make_grid_plots(dirname)
    elif cmd == "writeSteadyState":
        write_steady_state(dirname)
    elif cmd == "makeMDSImages":
        make_mds_images(dirname)
