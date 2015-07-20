#!/usr/bin/env python


"""

James McDermott [jmmcd@jmmcd.net]

This script runs some simple experiments with GA, GP and TSP operators.

If something doesn't work, please email me [jmmcd@jmmcd.net] and I
will help.
"""


import math
import numpy as np
import scipy.stats
import random
import sys
import os
import os.path
import itertools
from collections import OrderedDict
import matplotlib.pyplot as plt
import cPickle as pickle
from sklearn.manifold import MDS, TSNE

import random_walks
import plotting
import tsp
import bitstring

def operator_difference_and_compound_experiment():
    """Take the 6 permutation operators in pairs, and for each calculate
    the GINI/stddev/etc of the compound operator, and also the
    difference between the pair.

    """
    
    opss = ("two_opt", "twoh_opt", "three_opt", "three_opt_broad", "swap_two", "swap_adj")
    for n in range(6, 11):
        diff = np.zeros((len(opss), len(opss)))
        Gini = np.zeros((len(opss), len(opss)))
        SD = np.zeros((len(opss), len(opss)))
        CV = np.zeros((len(opss), len(opss)))
        for i, opi in enumerate(opss):
            for j, opj in enumerate(opss):
                basedir = sys.argv[1] + "/tsp_length_%d_" % n
                ops = [opi, opj]
                ps = [np.genfromtxt(basedir + op + "/TP_row0.dat") for op in ops]
                names = "+".join(ops)
                delta = random_walks.operator_difference_RMSE(*ps)
                print "delta", names, delta
                diff[i][j] = delta

                wts = [1.0 / len(ops) for _ in ops]
                compound_tp = random_walks.compound_operator(wts, ps)
                gini = random_walks.gini_coeff(compound_tp)
                sd = np.std(compound_tp)
                cv = np.std(compound_tp) / np.mean(compound_tp)
                print "compound gini", names, gini
                Gini[i][j] = gini
                print "compound stddev", names, sd
                SD[i][j] = sd
                print "compound coefvar", names, cv
                CV[i][j] = cv
                
        basedir = sys.argv[1] + "/permutation_size_%d_compound_diff/" % n
        
        filename = os.path.join(basedir + "diff")
        np.savetxt(filename+".dat", diff)
        plot_grid(filename, ("Operator difference %d" % n), opss, diff)
        
        filename = os.path.join(basedir + "diff_mds")
        op_diff_mds(filename, ("Operator difference %d" % n), opss, diff)

        filename = os.path.join(basedir + "compound_gini")
        np.savetxt(filename+".dat", Gini)
        plot_grid(filename, ("Compound operator Gini %d" % n), opss, Gini)

        filename = os.path.join(basedir + "compound_sd")
        np.savetxt(filename+".dat", SD)
        plot_grid(filename, ("Compound operator SD %d" % n), opss, SD)

        filename = os.path.join(basedir + "compound_cv")
        np.savetxt(filename+".dat", CV)
        plot_grid(filename, ("Compound operator CV %d" % n), opss, CV)
        
def op_diff_mds(filename, title, names, x):
    """Given a 2d array of operator differences, use MDS to lay them out and make aplot"""
    
    mds = MDS(dissimilarity="precomputed")
    npos = mds.fit_transform(x)
    # Make a figure
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    plt.scatter(npos[:,0], npos[:,1], s=40)

    names = map(lambda s: s.replace("_", " "), names)
    ax.set_xlim(np.min(npos[:,0]) * 1.25, np.max(npos[:,0]) * 1.25)
    ax.set_ylim(np.min(npos[:,1]) * 1.25, np.max(npos[:,1]) * 1.25)
    for label, x, y in zip(names, npos[:, 0], npos[:, 1]):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-10, 10),
            textcoords = 'offset points', ha = 'left', va = 'bottom')
            #bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            #arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
    plt.savefig(filename+".pdf")
    plt.savefig(filename+".eps")
    fig.clear()
    

def plot_grid(filename, title, names, x):
    """Given a 2d array, make a plot representing values by colour patches"""
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    plt.imshow(x, interpolation="nearest", cmap=plt.get_cmap("hot"))
    ax.set_aspect('equal')
    names = map(lambda s: s.replace("_", " "), names)
    plt.xticks(range(len(names)), names, rotation='vertical')
    plt.yticks(range(len(names)), names)
    plt.subplots_adjust(bottom=0.25)
    # ax.set_xticklabels(('',) + names)
    # ax.set_yticklabels(('',) + names)
    cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    plt.savefig(filename+".pdf")
    plt.savefig(filename+".eps")
    fig.clear()
            
def combinations_var_len(x):
    """Combinations of all lengths: 'ABC' -> '', 'A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC'"""
    for i in range(len(x) + 1):
        for item in itertools.combinations(x, i):
            yield item
            
def walk_permutation(n, op, nsteps):
    """Random walk on a permutation space"""
    result = []
    t = list(range(n))
    result.append(tuple(t))
    for i in range(nsteps):
        t = op(t)
        result.append(tuple(t))
    return result

def walk_bitstring(n, op, nsteps):
    result = []
    t = [0] * n
    result.append(tuple(t))
    for i in range(nsteps):
        t = op(t)
        result.append(tuple(t))
    return result

def walk_tree(n, op, nsteps):
    # call into Java for this?
    raise NotImplementedError()

def rw_experiment_with_tp(tp):
    N = len(tp)
    walklen = int(math.ceil(math.sqrt(N)))
        
    results = []
    fitvals = np.zeros_like(tp[0]) # don't need them for this
    reps = 30
    for rep in range(reps):
        samples, fit_samples, best = random_walks.hillclimb(tp, fitvals, walklen, rw=True)
        prop_unique = float(len(set(samples))) / len(samples)
        results.append(prop_unique)
    return np.mean(results), np.std(results)
    
    
def rw_experiment_with_op(space, n, op):
    """Proportion of unique individuals in a random walk: run many walks
     and calculate the proportion of unique individuals
     encountered. use a rw of length equal to sqrt of size of space.
     expect many duplicates even with explorative operators, but many
     more with exploitative ones

    """
    if space == "tree":
        raise NotImplementedError

    if space == "permutation":
        N = tsp.count_permutations(n)
    elif space == "bitstring":
        N = 2**n
    walklen = int(math.ceil(math.sqrt(N)))
    print "walklen", walklen
    results = []
    reps = 30
    for rep in range(reps):
        if space == "permutation":
            samples = walk_permutation(n, op, walklen)
        elif space == "bitstring":
            samples = walk_bitstring(n, op, walklen)
        elif space == "tree":
            raise NotImplementedError()

        x = float(len(set(samples))) / len(samples)
        results.append(x)
    print "XXX", np.mean(results), np.std(results)
    return np.mean(results), np.std(results)
    

def write_tp_row0(space):
    """Write out just the first row of the transition matrix for several
    operators of the given space. The first row is enough since later
    rows are just a rotation of the first.

    """
    
    basename = sys.argv[1]
    if space == "permutation":
        ops = ("two_opt", "twoh_opt", "three_opt", "three_opt_broad", "swap_two", "swap_adj")
        sizes = (9, 10)
    elif space == "bitstring":
        ops = ("per_gene", "per_ind")
        pmuts = (0.0001, 0.0003, 0.0010, 0.0033, 0.01, 0.0333, 0.1, 0.3333)
        sizes = (10, 12)
    for op in ops:
        for size in sizes:
            print op, size
            if space == "permutation":
                filename = os.path.join(basename,
                                        "space_%s/size_%d_op_%s" % (space, size, op),
                                        "TP_row0.dat")
                x = tsp.get_tm_first_row(size, move=op)
                np.savetxt(filename, x)
            elif space == "bitstring":
                if op == "per_ind":
                    filename = os.path.join(basename,
                                            "space_%s/size_%d_op_%s" % (space, size, op),
                                            "TP_row0.dat")
                    x = bitstring.generate_bitstring_tm_row0(size, pmut=None)
                    np.savetxt(filename, x)
                elif op == "per_gene":
                    for pmut in pmuts:
                        op_ = "%s_%.4f" % (op, pmut)
                        filename = os.path.join(basename,
                                                "space_%s/size_%d_op_%s" % (space, size, op_),
                                                "TP_row0.dat")
                        x = bitstring.generate_bitstring_tm_row0(size, pmut=pmut)
                        np.savetxt(filename, x)
                else:
                    raise ValueError("Unexpected operation + op")

def basic_stats_and_plots(space, do_rw_experiment=False):
    """For several operators, read in their first rows and calculate the
    basic stats (Gini, SD, CV) and nneighbours, and make the basic
    plots of these.

    """

    basename = sys.argv[1]

    if space == "permutation":
        ops = ("two_opt", "twoh_opt", "three_opt", "three_opt_broad", "swap_two", "swap_adj")
        opfs = {
            "two_opt": tsp.two_opt,
            "twoh_opt": tsp.twoh_opt,
            "three_opt": tsp.three_opt,
            "three_opt_broad": tsp.three_opt_broad,
            "swap_two": tsp.swap_two,
            "swap_adj": tsp.swap_adj
            }
        sizes = range(9, 11)
        
    elif space == "bitstring":
        ops = ("per_gene_0.0001", "per_gene_0.0003",
               "per_gene_0.0010", "per_gene_0.0033",
               "per_gene_0.0100", "per_gene_0.0333",
               "per_gene_0.1000", "per_gene_0.3333",
               "per_ind")
        opfs = {
            "per_ind": bitstring.bitstring_per_ind_mutation,
            "per_gene_0.0001": bitstring.make_bitstring_per_gene_mutation(0.0001),
            "per_gene_0.0003": bitstring.make_bitstring_per_gene_mutation(0.0003),
            "per_gene_0.0010": bitstring.make_bitstring_per_gene_mutation(0.0010),
            "per_gene_0.0033": bitstring.make_bitstring_per_gene_mutation(0.0033),
            "per_gene_0.0100": bitstring.make_bitstring_per_gene_mutation(0.0100),
            "per_gene_0.0333": bitstring.make_bitstring_per_gene_mutation(0.0333),
            "per_gene_0.1000": bitstring.make_bitstring_per_gene_mutation(0.1000),
            "per_gene_0.3333": bitstring.make_bitstring_per_gene_mutation(0.3333),
            }
        sizes = (10, 12)
    elif space == "tree" and do_rw_experiment:
        basic_stats_and_plots_tree()
        return

    prop_unique_v_sd_cv_gini = []
    for size in sizes:
        stddev = []
        coefvar = []
        gini = []
        nneighbours = []
        prop_unique = []
        for op in ops:
            filename = os.path.join(basename,
                                    "space_%s/size_%d_op_%s" % (space, size, op),
                                    "TP_row0.dat")
            print op, size
            x = np.genfromtxt(filename)
            # stats to get:
            sd = np.std(x)
            stddev.append(sd)
            cv = np.std(x)/np.mean(x)
            coefvar.append(cv)
            g = random_walks.gini_coeff(x)
            gini.append(g)
            nneighbours.append(np.sum(x > 0))
            if do_rw_experiment:
                mu, sigma = rw_experiment_with_op(space, size, opfs[op])
                prop_unique.append((mu, sigma))
                prop_unique_v_sd_cv_gini.append("%s %d %s %f %f %f %f %f" % (
                    space, size, op, sd, cv, g, mu, sigma))

        basename = sys.argv[1]

        # barcharts
        barchart(basename, space, size, "Gini", gini, ops)
        barchart(basename, space, size, "SD", stddev, ops)
        barchart(basename, space, size, "CV", coefvar, ops)

        # scatterplots v nneighbours
        scatterplot(basename, space, size, "Gini", "# Neighbours", gini, nneighbours, ops)
        scatterplot(basename, space, size, "SD", "# Neighbours", stddev, nneighbours, ops)
        scatterplot(basename, space, size, "CV", "# Neighbours", coefvar, nneighbours, ops)
        
        if do_rw_experiment:
            # scatterplots v prop unique
            scatterplot(basename, space, size, "Gini", "Prop unique", gini, prop_unique, ops)
            scatterplot(basename, space, size, "SD", "Prop unique", stddev, prop_unique, ops)
            scatterplot(basename, space, size, "CV", "Prop unique", coefvar, prop_unique, ops)

    if do_rw_experiment:
        filename = os.path.join(basename, "space_%s/prop_unique_v_sd_cv_gini.dat" % space)
        open(filename, "w").write("\n".join(prop_unique_v_sd_cv_gini))

def basic_stats_and_plots_tree():
    basename = sys.argv[1]
    filename = os.path.join(basename, "depth_2/TP.dat")
    tp = np.genfromtxt(filename)

    sd = random_walks.mu_sigma(tp)[0]
    g = random_walks.mean_gini_coeff(tp)
    cv = random_walks.mu_sigma_cv(tp)[0]
    mu, sigma = rw_experiment_with_tp(tp)
    
    prop_unique_v_sd_cv_gini = []
    space = "tree"
    op = "subtree"
    size = 2
    
    prop_unique_v_sd_cv_gini.append("%s %d %s %f %f %f %f %f" % (
        space, size, op, sd, cv, g, mu, sigma))
    filename = os.path.join(basename, "space_%s/prop_unique_v_sd_cv_gini.dat" % space)
    open(filename, "w").write("\n".join(prop_unique_v_sd_cv_gini))
                              
    

def scatterplot(basename, space, size, xlabel, ylabel, xs, ys, names):
    import inspect
    print "XXX", [locals()[arg] for arg in inspect.getargspec(scatterplot).args]
    
    names = map(lambda s: s.replace("_", " "), names)
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ylabel == "Prop unique":
        ys = np.array(ys)
        #plt.scatter(xs, ys[:,0], s=40)
        ax.errorbar(xs, ys[:,0], yerr=ys[:,1], fmt='o')
        annotate_ys = ys[:,1]
    else:
        plt.scatter(xs, ys, s=40)
        annotate_ys = ys
    # ax.set_xlim((0.999, 1.0002))
    for label, x, y in zip(names, xs, annotate_ys):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-10, 0),
            textcoords = 'offset points', ha = 'right', va = 'center')

    xlabel_s = xlabel.lower()
    if "Prop unique" in ylabel:
        ylabel_s = "prop_unique"
    else:
        ylabel_s = "nneighbours"
    filename = "%s/space_%s/size_%d_%s_v_%s" % (basename, space, size, xlabel_s, ylabel_s)
    plt.savefig(filename+".pdf")
    plt.savefig(filename+".eps")
    fig.clear()
    del ax
    del fig

def barchart(basename, space, size, ylabel, y, names):
    print "y", y
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)

    names = map(lambda s: s.replace("_", " "), names)
    ax.bar(range(len(names)), y, align='center', alpha=0.4)
    plt.xticks(range(len(names)), names, rotation='vertical')
    plt.ylabel(ylabel)
    plt.ylim((math.floor(100 * min(y)) / 100.0, math.ceil(100 * max(y)) / 100.0))
    #plt.ylim((0.998, 1.0))
    plt.subplots_adjust(bottom=0.4)
    ylabel_s = ylabel.lower()
    filename = "%s/space_%s/size_%d_%s_barchart" % (basename, space, size, ylabel_s)

    plt.savefig(filename+".pdf")
    plt.savefig(filename+".eps")
    fig.clear()
    del ax
    del fig

def plot_prop_unique_v_sd_cv_gini(spaces):
    plt.rcParams['legend.loc'] = 'best'
    for name, column in zip(["SD", "CV", "Gini"], [3, 4, 5]):
        basedir = sys.argv[1]
        fig = plt.figure(figsize=(6, 4.5))
        ax = fig.add_subplot(111)
        for space, sym, colour in zip(spaces, ("s", "o", "^"),
                              ["blue", "green", "red"]):
            x = np.genfromtxt(os.path.join(basedir,
                                           "space_"+space,
                                           "prop_unique_v_sd_cv_gini.dat"),
                              delimiter=" ").T
            plt.scatter(x[column], x[6], s=50, marker=sym, label=space, c=colour,
                        alpha=0.7)
        if name == "Gini":
            loc = 3 # lower left
        elif name == "CV":
            loc = 4 # lower right
        else:
            loc = 1
        plt.legend(loc=loc)
        plt.xlabel(name)
        plt.ylabel("Prop unique")
        filename = "%s/prop_unique_v_%s" % (basedir, name.lower())
        plt.savefig(filename+".pdf")
        plt.savefig(filename+".eps")
        fig.clear()
        del ax
        del fig
    
        
        

if __name__ == "__main__":
    spaces = ("bitstring", "permutation", "tree")
    # spaces = ("permutation", "bitstring")
    # spaces = ("tree",)
    for space in spaces:
        # if space != "tree":
        #     write_tp_row0(space)
        # basic_stats_and_plots(space, True)
        # if space == "permutation":
        #     operator_difference_and_compound_experiment()
        pass
    
    plot_prop_unique_v_sd_cv_gini(spaces)
