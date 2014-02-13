#!/usr/bin/env python


"""

James McDermott [jmmcd@jmmcd.net]

This script runs the simple experiments required for the EuroGP 2014
paper "Measuring mutation operators' exploration-exploitation
behaviour and long-term biases".

Get PODI from http://www.github.com/jmmcd/PODI.

Then adjust the paths in main() to suit your system.

Then run this function.

If something doesn't work, please email me [jmmcd@jmmcd.net] and I
will help.

Unfortunately it's all a bit hacky. Someday I will provide a simple
Makefile to regenerate each part, I hope.

"""


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

import random_walks
import plotting

def ga_hc_experiment(path_results):
    """Run some hill-climbs on variations of a GA space. Report
    performance."""
    uniformify_vals = [0.1, 0.5, .75, 0.9, 1.0, 1.0/0.9, 1.0/.75, 2.0, 10.0]
    noise_vals = [0, 1, 10, 100, 1000]
    results = OrderedDict()

    ga_length = 10
    tp_path = os.path.join(path_results, "ga_length_10", "TP.dat")
    try:
        ga_tp = np.genfromtxt(tp_path)
    except:
        ga_tp, _ = random_walks.generate_ga_tm(ga_length, pmut=1.0/ga_length)
        np.savetxt(tp_path, ga_tp)

    fit_path = os.path.join(path_results, "ga_length_10", "fitness_vals.dat")
    try:
        ga_fit = np.genfromtxt(fit_path)
    except:
        ga_fit = random_walks.onemax_fitvals(ga_length)
        np.savetxt(fit_path, ga_fit)

    # just get mu(sigma()), don't bother with sigma(sigma())
    mu_sigma_vals = [random_walks.mu_sigma(random_walks.uniformify(ga_tp, uniformify_val))[0]
                     for uniformify_val in uniformify_vals]

    reps = 30
    steps = 50
    for rep_name, tp, fitvals in [["ga", ga_tp, ga_fit]]:

        for noise_val in noise_vals:

            tmp_fit = random_walks.permute_vals(fitvals, noise_val)

            for uniformify_val in uniformify_vals:
                for rep in range(reps):
                    tp_tmp = random_walks.uniformify(tp, uniformify_val)
                    samples, fit_samples, best = random_walks.hillclimb(tp_tmp, tmp_fit,
                                                                        steps, rw=False)
                    x = best
                    results[rep_name, uniformify_val, noise_val, rep] = x
    return results, mu_sigma_vals

def plot_ga_hc_results(results, mu_sigma_vals, path_results):
    """Plot the results of the GA HC experiments above."""
    uniformify_vals = [0.1, 0.5, .75, 0.9, 1.0, 1.0/0.9, 1.0/.75, 2.0, 10.0]
    noise_vals = [0, 1, 10, 100, 1000]

    reps = 30
    for rep_name in ["ga"]:
        for noise_val in noise_vals:
            mu = []
            err = []

            for uniformify_val in uniformify_vals:
                x = [results[rep_name, uniformify_val, noise_val, rep]
                     for rep in range(reps)]
                mu.append(np.mean(x))
                err.append(np.std(x))

            plt.figure(figsize=(5, 2.5))
            plt.errorbar(mu_sigma_vals, mu, yerr=err, lw=3)
            plt.title(rep_name.upper() + r" OneMax with noise %d" % noise_val)
            plt.xlabel(r"$\mu(\sigma_r(p))$", fontsize=16)
            plt.ylabel("Fitness")
            plt.ylim(0, 10)
            filename = os.path.join(path_results, rep_name + "_noise_%d_hc" % noise_val)
            plt.savefig(filename + ".pdf")
            plt.savefig(filename + ".eps")

def write_ga_fitness_vals(n):
    ga_fit = onemax_fitvals(n)
    ga_outfile = open(dirname + "/ga_length_" + str(n) + "/fitness_vals.dat", "w")
    for fitval in ga_fit:
        ga_outfile.write(fitval)

def write_gp_trees(path_results):
    n = 2
    import generate_trees
    outfile = file(os.path.join(path_results, "depth_2", "all_trees.dat"), "w")
    trees_depths = generate_trees.trees_of_depth_LE(n,
                                                    ("x0", "x1"),
                                                    OrderedDict([("*", 2), ("+", 2),
                                                                 ("-", 2), ("AQ", 2)]),
                                                    as_string=False)
    for tree, depth in trees_depths:
        # hack hack: this is because a single variable gets a bare x0 instead of 'x0'.
        if len(tree) < 3:
            tree = '"%s"' % tree

        outfile.write(str(tree) + "\n")


def ga_gp_rw_experiment(path_results):
    uniformify_vals = [0.1, 0.5, .75, 0.9, 1.0, 1.0/0.9, 1.0/.75, 2.0, 10.0]
    results = OrderedDict()

    ga_length = 10
    gp_depth = 2

    tp_path = os.path.join(path_results, "ga_length_10", "TP.dat")
    ga_tp = np.genfromtxt(tp_path)

    tp_path = os.path.join(path_results, "depth_2", "TP.dat")
    gp_tp = np.genfromtxt(tp_path)

    print "doing probability of encounter experiment"
    # Do the "probability of encounter" experiment first
    reps = 100
    steps = 50

    fit_path = os.path.join(path_results, "ga_length_10", "fitness_vals.dat")
    ga_fit = np.genfromtxt(fit_path)

    gp_fit = [float(s) for s in
              open(os.path.join(path_results, "depth_2", "all_fitness_values.dat")).readlines()]

    inds = 0, len(gp_fit)-1
    hc_encounters = [0.0, 0.0]
    rw_encounters = [0.0, 0.0]
    for rep_name, tp, fitvals in [["gp", gp_tp, gp_fit]]:
        for rep in range(reps):
            samples, fit_samples, best = random_walks.hillclimb(gp_tp, fitvals, steps, rw=False)
            for i in range(2):
                if inds[i] in samples:
                    hc_encounters[i] += 1.0 / reps
            samples, fit_samples, best = random_walks.hillclimb(gp_tp, fitvals, steps, rw=True)
            for i in range(2):
                if inds[i] in samples:
                    rw_encounters[i] += 1.0 / reps
    print "hc_encounters", hc_encounters
    print "rw_encounters", rw_encounters

    # now the GA v GP hillclimb experiments
    reps = 30
    for rep_name, tp, fitvals in [["ga", ga_tp, ga_fit],
                                  ["gp", gp_tp, gp_fit]]:
        for uniformify_val in uniformify_vals:
            for rep in range(reps):
                tp_tmp = random_walks.uniformify(tp, uniformify_val)
                samples, fit_samples, best = random_walks.hillclimb(tp_tmp, fitvals, steps, rw=True)
                x = float(len(set(samples))) / len(samples)
                results[rep_name, uniformify_val, rep] = x
    return results, ga_fit, gp_fit


def plot_ga_gp_rw_results(results, mu_sigma_vals, path_results):
    uniformify_vals = [0.1, 0.5, .75, 0.9, 1.0, 1.0/0.9, 1.0/.75, 2.0, 10.0]

    reps = 30
    for rep_name in "ga", "gp":

        mu = []
        err = []
        for uniformify_val in uniformify_vals:
            x = [results[rep_name, uniformify_val, rep]
                 for rep in range(reps)]
            mu.append(np.mean(x))
            err.append(np.std(x))

        plt.figure(figsize=(5, 2.5))
        plt.errorbar(mu_sigma_vals, mu, yerr=err, lw=3)
        plt.title(rep_name.upper())
        plt.xlabel(r"$\mu(\sigma_r(p))$", fontsize=16)
        plt.ylabel("Exploration")
        plt.ylim(0, 1)
        filename = os.path.join(path_results, rep_name + "_uniformify_rw")
        plt.savefig(filename + ".pdf")
        plt.savefig(filename + ".eps")


def main()

    path_PODI_gp = "/tmp/PODI/src/gp.py"
    path_results = "/tmp/results/"

    # make results dir
    try:
        os.makedirs(os.path.join(path_results, "depth_2"))
    except OSError:
        pass
    try:
        os.makedirs(os.path.join(path_results, "ga_length_10"))
    except OSError:
        pass

    # compile and run the Java code for generating transition matrix
    cwd = os.getcwd()
    os.chdir("../../java")
    cmd = "make all"
    os.system(cmd)
    cmd = "make completeMatricesDepth2"
    os.system(cmd)
    os.chdir(cwd)

    generate all GP trees
    write_gp_trees(path_results)

    # use PODI's GP-fitness code to evaluate fitness of all GP trees
    cwd = os.getcwd()
    os.chdir(os.path.dirname(path_PODI_gp))
    sys.path.append(os.getcwd())
    import gp
    gp.read_trees_write_fitness_EuroGP2014(os.path.join(path_results, "depth_2", "all_trees.dat"),
                                           os.path.join(path_results, "depth_2", "all_fitness_values.dat"))
    os.chdir(cwd)

    results_file = os.path.join(path_results, "EuroGP_2014_results.pkl")
    try:
        # restore from a save, if it's been saved
        results = pickle.load(file(results_file))
        (ga_hc_results, ga_gp_rw_results, ga_fit, gp_fit, mu_sigma_vals) = results
    except:
        # run and save results
        ga_hc_results, mu_sigma_vals = ga_hc_experiment(path_results)
        ga_gp_rw_results, ga_fit, gp_fit = ga_gp_rw_experiment(path_results)
        results = (ga_hc_results, ga_gp_rw_results, ga_fit, gp_fit, mu_sigma_vals)
        pickle.dump(results, file(results_file, "w"))

    # plot
    plot_ga_hc_results(ga_hc_results, mu_sigma_vals, path_results)
    plot_ga_gp_rw_results(ga_gp_rw_results, mu_sigma_vals, path_results)


    # the stationary state and in-degree
    plotting.write_steady_state(os.path.join(path_results, "depth_2"))

if __name__ == "__main__":
    main()
