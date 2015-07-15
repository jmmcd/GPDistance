#!/usr/bin/env python


"""

James McDermott [jmmcd@jmmcd.net]

This script runs the simple experiments required for the EURO 2015
paper "Measuring exploration-exploitation behaviour of neighbourhood
operators in permutation search spaces".

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

def operator_difference_and_compound_experiment():
    """Take the 6 operators in pairs, and for each calculate the
    GINI/stddev of the compound operator, and also the difference
    between the pair.

    """
    
    opss = ("two_opt", "twoh_opt", "three_opt", "three_opt_broad", "swap", "swap_adj")
    for n in range(6, 11):
        diff_result = np.zeros((len(opss), len(opss)))
        comp_result = np.zeros((len(opss), len(opss)))
        comp_result_stddev = np.zeros((len(opss), len(opss)))
        for i, opi in enumerate(opss):
            for j, opj in enumerate(opss):
                basedir = "/Users/jmmcd/Dropbox/GPDistance/results/tsp_length_%d_" % n
                ops = [opi, opj]
                ps = [np.genfromtxt(basedir + op + "/TP_row0.dat") for op in ops]
                names = "+".join(ops)
                delta = random_walks.operator_difference_RMSE(*ps)
                print "delta", names, delta
                diff_result[i][j] = delta

                wts = [1.0 / len(ops) for _ in ops]
                compound_tp = random_walks.compound_operator(wts, ps)
                gini = random_walks.gini_coeff(compound_tp)
                stddev = np.std(compound_tp)
                print "compound gini", names, gini
                comp_result[i][j] = gini
                print "compound stddev", names, stddev
                comp_result_stddev[i][j] = stddev
                
        op_diff_filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_op_diff" % n
        np.savetxt(op_diff_filename+".dat", diff_result)
        plot_grid(op_diff_filename, ("Operator difference %d" % n), opss, diff_result)
        
        comp_op_filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_comp_op" % n
        np.savetxt(comp_op_filename+".dat", comp_result)
        plot_grid(comp_op_filename, ("Compound operator Gini %d" % n), opss, comp_result)

        comp_op_stddev_filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_comp_op_stddev" % n
        np.savetxt(comp_op_stddev_filename+".dat", comp_result_stddev)
        plot_grid(comp_op_stddev_filename, ("Compound operator stddev %d" % n), opss, comp_result_stddev)

        mds_filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_mds" % n
        op_diff_mds(mds_filename, ("Operator difference %d" % n), opss, diff_result)

def op_diff_mds(filename, title, names, x):
    """Given a 2d array of operator differences, use MDS to lay them out and make aplot"""
    
    mds = MDS(dissimilarity="precomputed")
    npos = mds.fit_transform(x)
    # Make a figure
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    plt.scatter(npos[:,0], npos[:,1], s=40)
    
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
    fig.clear()
    

def plot_grid(filename, title, names, x):
    """Given a 2d array, make a plot representing values by colour patches"""
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    plt.imshow(x, interpolation="nearest", cmap=plt.get_cmap("hot"))
    ax.set_aspect('equal')
    plt.xticks(range(6), names, rotation='vertical')
    plt.yticks(range(6), names)
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
    fig.clear()
            
def combinations_var_len(x):
    """Combinations of all lengths: 'ABC' -> '', 'A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC'"""
    for i in range(len(x) + 1):
        for item in itertools.combinations(x, i):
            yield item
            
def tsp_walk(n, op, nsteps):
    """Walk on a permutation space"""
    result = []
    t = list(range(n))
    result.append(tuple(t))
    for i in range(nsteps):
        t = op(t)
        result.append(tuple(t))
    return result
    
def rw_experiment_with_op(n, op):
    """Proportion of unique individuals in a random walk: run many walks
     and calculate the proportion of unique individuals
     encountered. use a rw of length equal to sqrt of size of space.
     expect many duplicates even with explorative operators, but many
     more with exploitative ones

    """
    N = tsp.count_permutations(n)
    walklen = int(math.ceil(math.sqrt(N)))
    print "walklen", walklen
    results = []
    reps = 30
    for rep in range(reps):
        samples = tsp_walk(n, op, walklen)
        x = float(len(set(samples))) / len(samples)
        results.append(x)
    return np.mean(results), np.std(results)
    

def write_tp_rows():
    """Write out just the first row of the transition matrix for several
    TSP operators. The first row is enough since later rows are just a
    rotation of the first.

    """
    
    basename = sys.argv[1]
    ops = ("two_opt", "twoh_opt", "three_opt", "three_opt_broad", "swap", "swap_adj")
    lengths = (6, 7, 8, 9, 10)
    for op in ops:
        for length in lengths:
            filename = os.path.join(basename,
                                    "tsp_length_%d_%s" % (length, op),
                                    "TP_row0.dat")
            print op, length
            x = tsp.get_tm_first_row(length, move=op)
            np.savetxt(filename, x)

def basic_stats_and_plots():
    """For several TSP operators, read in their first rows and calculate
    the basic stats (Gini and stddev) and nneighbours, and make the basic
    plots of these."""
    
    basename = sys.argv[1]
    ops = ("two_opt", "twoh_opt", "three_opt", "three_opt_broad", "swap", "swap_adj")
    opfs = {
        "two_opt": tsp.two_opt,
        "twoh_opt": tsp.twoh_opt,
        "three_opt": tsp.three_opt,
        "three_opt_broad": tsp.three_opt_broad,
        "swap": tsp.swap_two,
        "swap_adj": tsp.swap_adj
        }
    
    lengths = range(6, 11)
    for length in lengths:
        stddev = []
        gini = []
        nneighbours = []
        prop_unique = []
        for op in ops:
            filename = os.path.join(basename,
                                    "tsp_length_%d_%s" % (length, op),
                                    "TP_row0.dat")
            print op, length
            x = np.genfromtxt(filename)
            # stats to get:
            stddev.append(np.std(x))
            gini.append(random_walks.gini_coeff(x))
            nneighbours.append(np.sum(x > 0))
            mu, sigma = rw_experiment_with_op(length, opfs[op])
            prop_unique.append((mu, sigma))

        gini_barchart(length, gini, ops)
        stddev_barchart(length, stddev, ops)
        plot_gini_v_nneighbours(length, gini, nneighbours, ops)
        plot_stddev_v_nneighbours(length, stddev, nneighbours, ops)
        plot_gini_v_prop_unique(length, gini, prop_unique, ops)
        plot_stddev_v_prop_unique(length, stddev, prop_unique, ops)
        
def plot_gini_v_nneighbours(length, gini, nneighbours, ops):
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    ax.set_xlabel("Gini")
    ax.set_ylabel("# Neighbours")
    plt.scatter(gini, nneighbours, s=40)
    ax.set_xlim((0.999, 1.0002))
    for label, x, y in zip(ops, gini, nneighbours):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-10, 0),
            textcoords = 'offset points', ha = 'right', va = 'center')
        
    filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_gini_nneighbours" % length
    plt.savefig(filename+".pdf")
    fig.clear()

def plot_stddev_v_nneighbours(length, stddev, nneighbours, ops):
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    ax.set_xlabel("stddev")
    ax.set_ylabel("# Neighbours")
    plt.scatter(stddev, nneighbours, s=40)
    for label, x, y in zip(ops, stddev, nneighbours):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-10, 0),
            textcoords = 'offset points', ha = 'right', va = 'center')
    
    filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_stddev_nneighbours" % length
    plt.savefig(filename+".pdf")
    fig.clear()

def plot_gini_v_prop_unique(length, gini, prop_unique, ops):
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    ax.set_xlabel("Gini")
    ax.set_ylabel("Prop unique")
    prop_unique = np.array(prop_unique)
    ax.errorbar(gini, prop_unique[:,0], yerr=prop_unique[:,1], fmt='o')
    for label, x, y in zip(ops, gini, prop_unique[:, 0]):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-10, -10),
            textcoords = 'offset points', ha = 'right', va = 'center')
    filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_gini_prop_unique" % length
    # plt.subplots_adjust(bottom=0.15)
    plt.savefig(filename+".pdf")
    fig.clear()

def plot_stddev_v_prop_unique(length, stddev, prop_unique, ops):
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    ax.set_xlabel("stddev")
    ax.set_ylabel("Prop unique")
    prop_unique = np.array(prop_unique)
    ax.errorbar(stddev, prop_unique[:,0], yerr=prop_unique[:,1], fmt='o')
    for label, x, y in zip(ops, stddev, prop_unique[:, 0]):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-10, -10),
            textcoords = 'offset points', ha = 'right', va = 'center')

    filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_stddev_prop_unique" % length
    # plt.subplots_adjust(bottom=0.15)
    plt.savefig(filename+".pdf")
    fig.clear()

def gini_barchart(length, gini, names):
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    ax.bar(range(6), gini, align='center', alpha=0.4)
    plt.xticks(range(6), names, rotation='vertical')
    plt.ylabel('Gini')
    #plt.ylim((math.floor(100 * min(gini)) / 100.0, math.ceil(100 * max(gini)) / 100.0))
    plt.ylim((0.998, 1.0))
    plt.subplots_adjust(bottom=0.25)
    filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_gini_barchart" % length
    plt.savefig(filename+".pdf")
    fig.clear()

def stddev_barchart(length, stddev, names):
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    ax.bar(range(6), stddev, align='center', alpha=0.4)
    plt.xticks(range(6), names, rotation='vertical')
    plt.ylabel('stddev')
    plt.ylim((math.floor(100 * min(stddev)) / 100.0, math.ceil(100 * max(stddev)) / 100.0))
    plt.subplots_adjust(bottom=0.3)
    filename = "/Users/jmmcd/Dropbox/GPDistance/results/EURO2015/tsp_length_%d_stddev_barchart" % length
    plt.savefig(filename+".pdf")
    fig.clear()
    

if __name__ == "__main__":
    write_tp_rows()
    basic_stats_and_plots()
    operator_difference_and_compound_experiment()
