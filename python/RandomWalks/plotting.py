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
import scipy.stats
import scipy.stats.mstats
from random_walks import set_self_transition_zero

# MAXTICKS is 1000 in IndexLocator
class MyLocator(mpl.ticker.IndexLocator):
    MAXTICKS=1500

def graph_distance_names(dirname):
    if "depth_6" in dirname:
        return ["D_TP", "MFPT"], ["D$_\mathrm{TP}$", "MFPT"]
    else:
        return [
            "D_TP", "SD_TP",
            "MFPT", "CT", "MFPT_VLA", "CT_VLA",
            "RSP", "FE", "CT_amp",
            "SP", "STEPS",
            "MSTP_10"
            ], [
            "$D_{\mathrm{TP}}$", r"SD$_{\mathrm{TP}}$",
            "MFPT", "CT", "MFPT-VLA", "CT-VLA",
            "RSP", "FE", "CT-amp",
            "SP", "STEPS",
            "MSTP$_{10}$"
            ]

def syntactic_distance_names(dirname):
    if "ga_length" in dirname:
        return ["Hamming"]
    else:
        return [
            # "NodeCount", "MinDepth", "MeanDepth", "MaxDepth",
            # "Symmetry", "MeanFanout", "DiscreteMetric",
            "TED",
            "TAD0", "TAD1", "TAD2", "TAD3", "TAD4", "TAD5",
            "FVD",
            "NCD", 
            "OVD",
            ]

def make_grid_plots(dirname):
    if "depth" in dirname:
        if "depth_6" not in dirname:
            names = open(dirname + "/all_trees.dat").read().strip().split("\n")
            names = map(lambda x: x.strip(), names)
        else:
            names = None
    else:
        # Assume GA
        length = int(dirname.strip("/").split("_")[2])
        names = [bin(i)[2:] for i in range(2**length)]

    syn_names = syntactic_distance_names(dirname)
    grph_names, grph_tex_names = graph_distance_names(dirname)
        
    for name in grph_names + syn_names:
        w = np.genfromtxt(dirname + "/" + name + ".dat")
        if "depth_6" not in dirname:
            assert(len(w) == len(names))
        make_grid(w, False, dirname + "/" + name)
    print names # better to print them in a list somewhere than in the graph
        
def make_grid(w, names, filename):
    # we dont rescale the data. matshow() internally scales the data
    # so that the smallest numbers go to black and largest to white.

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


def metric_distortion(a, b):
    """Calculate the metric distortion between two metrics on the same
    space. a and b are samples of the distances between pairs of
    points (can be between all pairs).

    See Gupta's paper on embeddings: given a mapping f: X->Y (both
    metric spaces), the contraction of f is sup(d_X(x, y)/d_Y(f(x),
    f(y))); the expansion of f is sup(d_Y(f(x), f(y))/d_X(x, y)). The
    distortion is contraction * expansion. The best (lowest) possible
    value is 1. To avoid divide-by-zero, this assumes d(x, y) != 0 --
    ie x != y, and d is strictly positive for distinct elements, which
    is true of any metric but not true of some of our distance
    functions. Below, we overwrite any divide-by-zeros just in case.

    In our case, X = Y = the search space, eg trees of depth 2 or
    less. d_X is one metric, eg TED; d_Y is another, eg D_TP; f is the
    identity mapping. This leads to these equations:

    contraction = sup(TED(x, y)/TP(x, y)); expansion = sup(TP(x,
    y)/TED(x, y)).

    Also, it means that the order of arguments (a, b) is
    unimportant."""

    old = np.seterr(divide="ignore")
    ab = a/b
    ab[~np.isfinite(ab)] = -1
    contraction = np.max(ab)

    ba = b/a
    ba[~np.isfinite(ba)] = -1
    expansion = np.max(ba)
    np.seterr(**old)
    return expansion * contraction

def metric_distortion_agreement(a, b):
    # because we want a measure of agreement
    return 1.0 / metric_distortion(a, b)
    
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

    # make sure we raise any error (so we can catch it), don't just
    # splat it on the terminal
    old = np.seterr(all='raise')
    try:
        if (isinstance(x, np.ma.core.MaskedArray) or
            isinstance(y, np.ma.core.MaskedArray)):
            corr, p = scipy.stats.mstats.kendalltau(x, y)
        else:
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

    # Make sure we raise any error (so we can catch it), don't just
    # splat it on the terminal. However we will ignore underflow
    # because it seems to happen with every single one. Similar
    # approach here:
    # [https://code.google.com/p/wnd-charm/source/browse/pychrm/trunk/pychrm/FeatureSet.py?r=723]
    old = np.seterr(all='raise', under='ignore')
    try:
        if (isinstance(x, np.ma.core.MaskedArray) or
            isinstance(y, np.ma.core.MaskedArray)):
            corr, p = scipy.stats.mstats.spearmanr(x, y)
        else:
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

    # Make sure we raise any error (so we can catch it), don't just
    # splat it on the terminal. 
    old = np.seterr(all='raise')
    try:
        if (isinstance(x, np.ma.core.MaskedArray) or
            isinstance(y, np.ma.core.MaskedArray)):
            corr, p = scipy.stats.mstats.pearsonr(x, y)
        else:
            corr, p = scipy.stats.pearsonr(x, y)
    except FloatingPointError:
        corr, p = 0.0, 1.0
    # restore old error settings
    np.seterr(**old)
    return corr, p

def make_correlation_tables(dirname, txt=""):

    syn_names = syntactic_distance_names(dirname)
    grph_names, grph_tex_names = graph_distance_names(dirname)
        
    d = load_data_and_reshape(dirname, syn_names + grph_names)

    def do_line(dist, dist_name):
        line = dist_name.replace("_TP", r"$_{\mathrm{TP}}$")
        for graph_distance in grph_names:
            print("getting association between " + graph_distance + " " + dist)
            if corr_type == "spearmanrho":
                corr, p = get_spearman_rho(d[graph_distance], d[dist])
            elif corr_type == "kendalltau":
                corr, p = get_kendall_tau(d[graph_distance], d[dist])
            elif corr_type == "pearsonr":
                corr, p = get_pearson_r(d[graph_distance], d[dist])
            elif corr_type == "metric_distortion":
                # there is no p-value associated with metric
                # distortion so use 1.0
                corr, p = metric_distortion_agreement(d[graph_distance], d[dist]), 1.0
            else:
                print("Unknown correlation type " + corr_type)
            if corr > 0.0 and p < 0.05:
                sig = "*"
            else:
                sig = " "
            line += r" & {0:1.2f} \hfill {1} ".format(corr, sig)
        line += r"\\"
        f.write(line + "\n")
        
    for corr_type in ["spearmanrho", "pearsonr", "kendalltau", "metric_distortion"]:
        if corr_type == "kendalltau" and len(d["D_TP"]) > 1000:
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
        elif corr_type == "metric_distortion":
            f.write("using metric distortion: ")
            
        f.write(txt)
        f.write(r"""\label{tab:correlationresults_"""
                + os.path.basename(dirname)
                + "_" + corr_type + "}}\n")
        f.write(r"""\begin{tabular}{""" + "|".join("l" for i in range(1+len(grph_names))) + r"""}
     & """ + " & ".join(grph_tex_names)  + r"""\\
\hline
\hline
""")

        for name, tex_name in zip(grph_names, grph_tex_names):
            do_line(name, tex_name)
        f.write(r"""\hline
\hline
""")
        for syn in syn_names:
            do_line(syn, syn)

        f.write(r"""\end{tabular}
\end{table*}
""")
        f.close()

        
def load_data_and_reshape(dirname, names):
    d = {}
    for name in names:
        # print("reading " + name)
        if "estimate_MFPT" in dirname:
            m = np.genfromtxt(dirname + "/" + name + ".dat", usemask=True, missing_values="NaN")
        else:
            m = np.genfromtxt(dirname + "/" + name + ".dat")
        d[name] = m.reshape(len(m)**2)
    return d

def make_scatter_plots(dirname):
    syn_names = syntactic_distance_names(dirname)
    grph_names, grph_tex_names = graph_distance_names(dirname)
        
    d = load_data_and_reshape(dirname, syn_names + grph_names)

    # while we have the data loaded in d, make scatter plots

    # graph v graph first
    for name1, tex_name1 in zip(grph_names, grph_tex_names):
        for name2, tex_name2 in zip(grph_names, grph_tex_names):
            if name1 < name2:
                # avoid plotting anything against itself, or plotting any pair twice
                make_scatter_plot(dirname, d, name1, tex_name1, name2, tex_name2)
                
    # graph v syn
    for name1, tex_name1 in zip(syn_names, syn_names):
        for name2, tex_name2 in zip(grph_names, grph_tex_names):
            make_scatter_plot(dirname, d, name1, tex_name1, name2, tex_name2)

def make_scatter_plot(dirname, d, name1, tex_name1, name2, tex_name2):
    filename = dirname + "/scatter_" + name1 + "_" + name2
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(d[name1], d[name2])
    ax.set_xlabel(tex_name1)
    ax.set_ylabel(tex_name2)
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png")

def make_histograms(dirname):
    syn_names = syntactic_distance_names(dirname)
    grph_names, grph_tex_names = graph_distance_names(dirname)
        
    d = load_data_and_reshape(dirname, syn_names + grph_names)

    # graph and syn names
    for name, tex_name in zip(grph_names + syn_names, grph_tex_names + syn_names):
        make_histogram(dirname, d, name, tex_name)
                
def make_histogram(dirname, d, name, tex_name, marks):
    infilename = dirname + "/" + name + ".dat"
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(d[name], 20, label=tex_name)
    ax.legend()
    fig.savefig(dirname + "/histogram_" + name + ".pdf")

def compare_TP_estimate_v_exact(dirname):
    
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

def compare_MFPT_estimate_RW_v_exact(dirname):

    def get_indices_of_common_entries(a, b):
        result = []
        for i in b:
            try:
                j = a.index(i)
                result.append(j)
            except ValueError:
                pass
        assert len(result) == len(b)
        return result
    
    filename = dirname + "/MFPT.dat"
    mfpt = np.genfromtxt(filename)
        
    filename = dirname + "/compare_MFPT_estimate_RW_v_exact.tex"
    f = open(filename, "w")

    if "depth_2" in dirname:
        lengths = [1298, 12980, 129800]
    elif "depth_1" in dirname:
        # NB with 18, too few values to run correlation
        lengths = [180, 1800, 18000]
    for length in lengths:
        
        # mfpte: read, mask nan, mask len < 5 (100x100)
        filename = dirname + "/estimate_MFPT_using_RW_" + str(length) + "/MFPT.dat"
        mfpte = np.genfromtxt(filename, usemask=True, missing_values="NaN")
        filename = dirname + "/estimate_MFPT_using_RW_" + str(length) + "/MFPT_len.dat"
        mfpte_len = np.genfromtxt(filename)
        min_vals = 5 # an attempt at reliability
        print("%d of %d values of length < 5" % (np.sum(mfpte_len < min_vals), len(mfpte_len)**2))
        mfpte[mfpte_len < min_vals] = np.ma.masked

        # mfpt: copy, select sampled only to make it 100x100
        mfpt_tmp = mfpt.copy()
        # need to restrict mfpt_tmp to the 100x100 entries which are
        # indicated by the trees_sampled.dat file
        filename = dirname + "/estimate_MFPT_using_RW_" + str(length) + "/trees_sampled.dat"
        trees_sampled = open(filename).read().strip().split("\n")
        filename = dirname + "/all_trees.dat"
        all_trees = open(filename).read().strip().split("\n")
        indices = get_indices_of_common_entries(all_trees, trees_sampled)
        # the selected indices are into both the rows and columns
        mfpt_tmp = mfpt_tmp[indices][:,indices]

        # mfpte will contain the self-hitting time on the diagonal: we
        # want zero there for true comparison.
        set_self_transition_zero(mfpte)

        # reshape both
        mfpte = mfpte.reshape(len(mfpte)**2)
        mfpt_tmp = mfpt_tmp.reshape(len(mfpt_tmp)**2)

        # correlate using mask
        f.write("Number of samples: " + str(length) + "\n")
        corr, p = get_pearson_r(mfpt_tmp, mfpte)
        f.write("Pearson R correlation " + str(corr) + "; ")
        f.write("p-value " + str(p) + ". ")
        corr, p = get_spearman_rho(mfpt_tmp, mfpte)
        f.write("Spearman rho correlation " + str(corr) + "; ")
        f.write("p-value " + str(p) + ". ")
        if len(mfpte) < 1000:
            corr, p = get_kendall_tau(mfpt_tmp, mfpte)
            f.write("Kendall tau correlation " + str(corr) + "; ")
            f.write("p-value " + str(p) + ". ")
        else:
            f.write("Omitting Kendall tau because it is infeasible for large matrices. ")
        f.write("\n")
    f.close()
    
def write_steady_state(dirname):
    """Get steady state, write it out and plot it. Also calculate the
    in-degree of each node, by summing columns. Plot that on the same
    plot. Calculate the correlation between steady-state and
    in-degree. Calculate the stddev of the steady-state as well, and
    (why not) the stddev of the TP matrix as well."""
    from random_walks import get_steady_state
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

    if "depth" in dirname:
        names = ["CT", "SD_TP", "TED", "TAD0", "OVD", "FVD"]
        # read in the trees
        filename = dirname + "/all_trees.dat"
        labels = open(filename).read().strip().split("\n")
    else:
        names = ["CT", "SD_TP", "Hamming"]
        labels = None
    for name in names:
        m = np.genfromtxt(dirname + "/" + name + ".dat")
        make_mds_image(m, dirname + "/" + name + "_MDS", labels)
    
def make_mds_image(m, filename, labels=None):
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

    if labels != None:
        # hard-coded to depth-2
        indices = [0, 2, 50, 52]
        for i in indices:
            plt.text(p[i,0], p[i,1], labels[i], style='italic',
                    bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
    
    plt.savefig(filename + ".png")
    plt.savefig(filename + ".pdf")

    # Could use custom marker types and colours to get an interesting
    # diagram?
    # http://matplotlib.org/api/path_api.html#matplotlib.path.Path
    # (pass marker=Path() from above), and then use random colours?

if __name__ == "__main__":
    cmd = sys.argv[1]
    dirname = sys.argv[2]

    if cmd == "compareTPEstimateVExact":
        compare_TP_estimate_v_exact(dirname)
    elif cmd == "compareMFPTEstimateRWVExact":
        compare_MFPT_estimate_RW_v_exact(dirname)
    elif cmd == "makeCorrelationTable":
        txt = sys.argv[3]
        make_correlation_tables(dirname, txt)
    elif cmd == "makeGridPlots": 
        make_grid_plots(dirname)
    elif cmd == "writeSteadyState":
        write_steady_state(dirname)
    elif cmd == "makeMDSImages":
        make_mds_images(dirname)
    elif cmd == "makeScatterPlots":
        make_scatter_plots(dirname)
    elif cmd == "makeHistograms":
        make_histograms(dirname)
