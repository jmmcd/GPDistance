package GP;

import TreeDistance.*;

import java.util.*;
import java.io.*;

import static java.lang.Math.*;


/*
 * Sample.java. This is the top-level class for doing sampling of
 * individuals and calculating syntactic and mutation-based distances
 * between them.
 *
 */

public class Sample {

    Random rng;
    int maxDepth;
    Language language;
    Mutation mutator;

    double epsilon = 0.000001;
    
    public Sample(int _maxDepth) {

        language = new Language(_maxDepth);
        rng = new Random();
        maxDepth = _maxDepth;
        mutator = new Mutation(language, maxDepth, rng);

    }

    public double mean(ArrayList<Integer> a) {
        assert a.size() > 0;
        double sum = 0.0;
        for (Integer i: a) {
            sum += i;
        }
        return sum / a.size();
    }

    // write out a matrix of values
    public void writeMatrix(double [][]vals,
                            String dirname,
                            String filename
                            ) {
        try {
            (new File(dirname)).mkdirs();

            // Open file
            FileWriter fw = new FileWriter(dirname + "/" + filename);

            for (int s = 0; s < vals.length; s++) {
                for (int t = 0; t < vals[s].length; t++) {
                    fw.write(vals[s][t] + " ");
                }
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    // Given a list of trees, write out a matrix of distances between
    // all pairs, for many different types of distance.
    public void writeMatrices(ArrayList<String> trees,
                              String dirname,
                              boolean writeNonNormalised
                              ) {

        String distanceNames[] = {
            "TP", "NCD", "FVD",
            "NodeCount", "MinDepth", "MeanDepth", "MaxDepth",
            "Symmetry", "MeanFanout", "DiscreteMetric",
            "TED",
            "TAD0", "TAD1", "TAD2", "TAD3", "TAD4", "TAD5",
            "OVD"
        };

        int n = trees.size();
        
        HashMap<String, double[][]> valss = new HashMap<String, double[][]>();
        for (String distance: distanceNames) {
            valss.put(distance, new double[n][n]);
        }

        // for every source tree
        for (int si = 0; si < trees.size(); si++) {
            String s = trees.get(si);
            
            // for every destination tree
            for (int ti = 0; ti < trees.size(); ti++) {
                String t = trees.get(ti);
                 
                // Get distances...
                HashMap<String, Double> distances = getDistances(s, t);
                
                // put them into the matrices
                for (String distance: distanceNames) {
                    valss.get(distance)[si][ti] = distances.get(distance);
                }
            }
        }

        for (String distance: distanceNames) {

            String filename;
            if (distance.equals("TP") && writeNonNormalised) {
                // Since we're sampling from the space, the
                // transition probabilities won't sum to 1. So
                // save to _nonnormalised.
                filename = "TP_nonnormalised.dat";
            } else {
                filename = distance + ".dat";
            }
            
            writeMatrix(valss.get(distance), dirname, filename);
        }
            
    }

    // This is non-uniform sampling
    public ArrayList<String> sampleByGrow(int n) {

        ArrayList<String> retval = new ArrayList<String>();
        while (retval.size() < n) {
            Tree tree = new Tree("x");
            mutator.grow(tree.getRoot(), maxDepth);
            String s = tree.toString();
            if (retval.indexOf(s) == -1) {
                retval.add(s);
            }
        }
        return retval;
    }

    // For each individual, calculate many mutations, and make an
    // estimate of the transition probabilities
    public void sampleOneStepProbabilities(ArrayList<String> trees,
                                                        int n,
                                                        String filename) {
        try {
            FileWriter file = new FileWriter(filename);
            for (String srcS: trees) {
                System.out.println("Working on " + srcS + " (" + trees.indexOf(srcS) + " of " + trees.size() + ")");
                // make a zero array to hold values
                double[] counts = new double[trees.size()];
                for (int i = 0; i < counts.length; i++) {
                    counts[i] = 0;
                }
                Tree src = new Tree(srcS);
                // perform many mutations and see where they go
                for (int i = 0; i < n * trees.size(); i++) {
                    // mutate src
                    Tree dest = mutator.mutate(src);
                    String destS = dest.toString();
                    // find new ind in the list
                    int idx = trees.indexOf(destS);
                    // +1 to new index
                    counts[idx]++;
                }
                // write normalised counts
                double sum = 0.0;
                for (int i = 0; i < trees.size(); i++) {
                    sum += counts[i];
                }
                for (int i = 0; i < trees.size(); i++) {
                    file.write(counts[i] / sum + " ");
                }
                // write a newline
                file.write("\n");
            }

            // close file
            file.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return;
    }
    

    // carry out the estimation of transition probabilities multiple
    // times for multiple different samples.
    public void estimateTPWithSuperNode(String basename, int reps, int N) {
        (new File(basename)).mkdirs();
        for (int i = 0; i < reps; i++) {
            int nNeighbours = (N - 2) / 2;
            // estimateTPBetweenPairWithSuperNode(basename,
            //                                    samplePairAndNeighbours(nNeighbours));
            estimateTPBetweenPairWithSuperNode(basename,
                                               sampleRandomWalk(5, 100, 10));
        }
    }

    // Estimate transition probabilities by taking a sample of
    // individuals and then modelling the rest of the space as a
    // single "super node". Can calculate the TP into the super node
    // as 1 minus the sum of all other TPs, etc. But we don't really
    // care about all these TPs, we really want to estimate MFPT. So
    // call out to the Python/Numpy code to estimate that. 
    public void estimateTPBetweenPairWithSuperNode(String basename,
                                                   ArrayList<String> sample) {

        int L = sample.size();
        
        try {

            FileWriter tpFile = new FileWriter(basename + "/TP.dat");
        
            double meanRowsum = 0.0; // sum of inward TPs to supernode over all rows but last
            for (int i = 0; i < L; i++) {
                double rowsum = 0.0;
                for (int j = 0; j < L; j++) {
                    double tp = mutator.transitionProbability(new Tree(sample.get(i)),
                                                              new Tree(sample.get(j)));
                    rowsum += tp;
                    tpFile.write(tp + " ");

                    // System.out.println("From " + sample.get(i));
                    // System.out.println("To " + sample.get(j));
                    // System.out.println("TP " + tp);

                }

                // write the inward TP for supernode
                tpFile.write(1.0 - rowsum + "");
                meanRowsum += (1.0 - rowsum);

                // System.out.println("rowsum " + rowsum);
                tpFile.write("\n");
            }

            meanRowsum /= L;

            // TP(x, S) = 1 - sum(TP(x, y)) for all y = 1 - rowsum
            
            // assume TP(S, S) = mean(rowsum)

            // TP(S, x) = (1 - mean(rowsum)) / L -- share out remaining probability

            for (int i = 0; i < L; i++) {
                tpFile.write((1.0 - meanRowsum) / L + " ");
            }
            tpFile.write(meanRowsum + "\n");

            // close file
            tpFile.close();
        
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Call Python MFPT code. It will read the TP matrix we just
        // wrote, then write out a new MFPT matrix to the same dir. 
        try {
            Process p = Runtime.getRuntime().exec("python ../python/RandomWalks/random_walks.py " + basename);
            p.waitFor();
            // BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            // String line=reader.readLine();

        } catch (Exception e) {
            System.out.println("Error " + e.getMessage());
        }

        // We can then read estimates of MFPT from it. We'll only read
        // the first row, and exclude node 0 -> node 0, and exclude
        // node 0 -> supernode.
        try {
            FileReader mfptFile = new FileReader(basename + "/MFPT.dat");
            LineNumberReader ln = new LineNumberReader(mfptFile);
            String line = ln.readLine();
            String[] vals = line.split(" ");

            for (int i = 1; i < vals.length - 1; i++) {
                String val = vals[i];
                HashMap<String, Double> distances = getDistances(sample.get(0),
                                                                 sample.get(i));
                System.out.println(val + " " + distances.get("TED") + " " + distances.get("FVD") + " " + distances.get("NCD"));
            }
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }

    }

    // Estimate the length of a random walk by simulation. This
    // function doesn't collect a sample of individuals. It just
    // performs walks and saves the number of steps between pairs of
    // individuals. So it aims to estimate FMPT by simulation rather
    // than exact methods on a sample from the space.

    // I wrote up a blog post explaining this algorithm here:
    // http://jmmcd.net/2013/09/28/simulating-random-walks.html

    // The basic idea is: Each time we hit an individual of interest
    // t, we can save a sampled walklength(s, t) = lastOccur(t) -
    // lastOccur(s) if we have previously seen s. The picture is:
    // others -> s -> others -> t. Note that it can happen that we see
    // s multiple times before seeing t: s0 -> others -> s1 -> others
    // -> s2 -> others -> t. In that case we should save (t-s0), but
    // we should not save (t-s1) or (t-s2). We can also see t multiple
    // times: s -> others -> t0 -> others -> t1 -> others -> t2. In
    // that case we save (t0-s), but we should not save (t1-s) or
    // (t2-s).

    // We can accomplish this by good book-keeping. They key idea is
    // we can't be in a walk from u to v and from v to u at the same
    // time. Can't be in two walks from u to v at the same time. the
    // variable rwStarted[u][v] means the time-step at which a walk
    // from u to v started. If we encounter u again before v, that is
    // *not* the start of a new walk from u to v. We save a sample
    // only when a walk ends on v, and then record that a new walk is
    // starting from v.
    
    public void randomWalking(
                              int nsteps,
                              ArrayList<String> selected,
                              int nsaves
                              ) {

        int n = selected.size();

        // for each pair of individuals (u, v) in selected, rwStarted
        // holds the time-step at which the rw from u to v began. if
        // we see i several times in a row before seeing j, then only
        // the *first* is relevant
        long[][] rwStarted = new long[n][n];
        int[][] placeToSave = new int[n][n];
        
        // for each pair of individuals in selected, walkLengths
        // holds samples of walkLengths between them.
        long[][][] walkLengths = new long[n][n][nsaves];

        // set up a hashmap for fast access
        HashMap<String, Integer> stringToIndex = new HashMap<String, Integer>();
        for (int i = 0; i < n; i++) {
            stringToIndex.put(selected.get(i), i);
        }

        // set up various matrices
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                rwStarted[i][j] = -1;
                placeToSave[i][j] = 0;
                for (int k = 0; k < nsaves; k++) {
                    walkLengths[i][j][k] = -1;
                }
            }
        }

        // start at one of the individuals of interest
        String sv = selected.get(0);
        Tree tv = new Tree(sv);

        for (long t = 0; t < nsteps; t += 1) {
            // System.out.println("");
            
            Integer v = stringToIndex.get(sv);
            if (v != null) {
                // we have found an individual sv of interest
                
                for (int u = 0; u < n; u++) {

                    if (rwStarted[u][v] == rwStarted[v][u] &&
                        rwStarted[u][v] == -1) {
                        // never saw u or v before, but this is the
                        // beginning of a walk from v to u
                        rwStarted[v][u] = t;

                    } else if (rwStarted[u][v] > -1) {
                        // we are currently in a walk from u to v,
                        // have now reached v, so if we haven't
                        // already collected plenty of (u, v) samples
                        // then collect this one
                        int i = placeToSave[u][v];
                        if (i < nsaves) {
                            walkLengths[u][v][i] = (t - rwStarted[u][v]);
                            placeToSave[u][v] += 1;
                        }

                        // now start a walk from v to u. note the
                        // order of these two assignments: does the
                        // right thing in the case of v == u, ie on
                        // diagonal
                        rwStarted[u][v] = -1;
                        rwStarted[v][u] = t;
                        
                    } else if (rwStarted[v][u] > -1) {
                        // we are in a walk from v to u, and have
                        // re-encountered v. ignore.
                        
                    } else {
                        System.out.println("Unexpected, in a walk from u and v and v to u at the same time?");
                        System.exit(1);
                    }
                }
            }

            tv = mutator.mutate(tv);
            sv = tv.toString();

            if (t % 10000 == 0) {
                System.out.println("" + t + " of " + nsteps + " mutations done");
            }
        }

        double mfpte[][] = new double[n][n];
        double mfpte_stddev[][] = new double[n][n];
        double mfpte_len[][] = new double[n][n];

        for (int uidx = 0; uidx < n; uidx++) {
            String u = selected.get(uidx);
            for (int vidx = 0; vidx < n; vidx++) {
                String v = selected.get(vidx);

                int nsaved = 0;
                double sum = 0.0;
                for (Long i: walkLengths[uidx][vidx]) {
                    if (i > -1) {
                        nsaved += 1;
                        sum += i;
                    }
                }
                double mean = sum / nsaved;
                double var = 0.0;
                for (Long i: walkLengths[uidx][vidx]) {
                    if (i > -1) {
                        var += (i - mean) * (i - mean);
                    }
                }
                var /= nsaved;
                var = sqrt(var);
                
                mfpte[uidx][vidx] = mean;
                mfpte_stddev[uidx][vidx] = var;
                mfpte_len[uidx][vidx] = (double) nsaved;
            }
        }

        String dirname = "../results/depth_" + maxDepth + "/estimate_MFPT_using_RW_" + nsteps;
        writeMatrices(selected,
                      dirname,
                      true);
        writeMatrix(mfpte, dirname, "/MFPTE.dat");
        writeMatrix(mfpte_stddev, dirname, "/MFPTE_stddev.dat");
        writeMatrix(mfpte_len, dirname, "/MFPTE_len.dat");

        writeListOfTrees(selected, dirname, "/trees_sampled.dat");
    }

    public void writeListOfTrees(ArrayList<String> trees,
                                 String dirname,
                                 String filename) {
        try {
            (new File(dirname)).mkdirs();

            // Open file
            FileWriter fw = new FileWriter(dirname + "/" + filename);

            for (String tree: trees) {
                fw.write(tree + "\n");
            }
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
            
            

    // Sample a pair and some neighbours. M is the number of
    // neighbours of each
    public ArrayList<String> samplePairAndNeighbours
        (int M) {

        // Get 2 trees and M neighbours of each into a list
        int N = 2;
        
        int L = N + N * M; // number of trees in the list
        
        ArrayList<Tree> ijAndNeighbours = new ArrayList<Tree>();
        ArrayList<String> ijAndNeighboursStrings = new ArrayList<String>();
        
        for (int i = 0; i < N; i++) {
            Tree tree = new Tree("x");
            mutator.grow(tree.getRoot(), maxDepth);
            ijAndNeighbours.add(tree.clone());
            ijAndNeighboursStrings.add(tree.toString());
        }
        
        for (int i = 0; i < N; i++) {
            while (ijAndNeighbours.size() < L) {
                Tree tree = ijAndNeighbours.get(i);
                // FIXME this is a horrible bug: have to clone() here,
                // else the TP can be zero or otherwise wrong. I think
                // because the TP algorithm gets confused about
                // ancestor of root.
                Tree newTree = mutator.mutate(tree).clone();
                if (ijAndNeighboursStrings.indexOf(newTree.toString()) == -1) {
                    ijAndNeighbours.add(newTree);
                    ijAndNeighboursStrings.add(newTree.toString());
                }
            }
        }

        return ijAndNeighboursStrings;
    }
        


    // Sample by random-walking. Aim is to collect a sample of size
    // precisely nindividuals. Use multiple walks to do so if needed,
    // each lasting up to nsteps. The distances from the first
    // individual to all the others will probably be the most
    // reliable, so when we get the MFPT matrix (for example), we
    // should only look at the top row.
    public ArrayList<String> sampleRandomWalk
        (int nsteps, int nindividuals, int ntries) {

        ArrayList<String> retval = new ArrayList<String>();

        Tree gamma_0 = new Tree("x");
        mutator.grow(gamma_0.getRoot(), maxDepth);
        retval.add(gamma_0.toString());

        while (retval.size() < nindividuals) {
            Tree gamma_i = gamma_0;

            for (int j = 0; j < nsteps && retval.size() < nindividuals; j++) {

                // Try up to ntries times to find an individual not
                // already in the sample. If we have to give up, don't
                // add anything on this iteration.
                for (int k = 0; k < ntries; k++) {
                    gamma_i = mutator.mutate(gamma_i);
                    String s = gamma_i.toString();
                    if (retval.indexOf(s) == -1) {
                        retval.add(s);
                        break;
                    }
                }
            }
        }
        return retval;
    }


    // This does Metropolis-Hastings sampling more or less as adapted
    // for GP by Vanneschi (see Vanneschi PhD thesis available online
    // [http://personal.disco.unimib.it/Vanneschi/thesis_vanneschi.pdf],
    // p 130). We run M-H multiple times if needed, starting from the
    // same "centre" node, until we have collected nindividuals, then
    // return the list. The distances from the first individual to all
    // the others will probably be the most reliable, so when we get
    // the MFPT matrix (for example), we should only look at the top
    // row.
    public ArrayList<String>
        sampleMetropolisHastings
        (int nsteps, int nindividuals, int ntries) {

        ArrayList<String> retval = new ArrayList<String>();

        Tree gamma_0 = new Tree("x");
        mutator.grow(gamma_0.getRoot(), maxDepth);
        float fgamma_0 = (float) (gamma_0.fitness());
        retval.add(gamma_0.toString());

        while (retval.size() < nindividuals) {
            Tree gamma_i = gamma_0;
            float fgamma_i = fgamma_0;

            for (int j = 0; j < nsteps && retval.size() < nindividuals; j++) {

                // Try up to ntries times to find a *better*
                // individual not already in the sample. 
                for (int k = 0; true; k++) {
                    Tree delta = mutator.mutate(gamma_i);
                    float fdelta = (float) delta.fitness();
                    double alpha;
                    // FIXME maximise should be passed-in as a
                    // parameter according to the problem.
                    boolean maximise = false;
                    if (maximise) {
                        alpha = min(1.0f, fdelta / fgamma_i);
                    } else {
                        alpha = min(1.0f, fgamma_i / fdelta);
                    }
                    if (Double.isNaN(alpha) || alpha <= epsilon) {
                        // It can happen. It's not worth worrying about.
                        alpha = epsilon;
                    }
                    // System.out.println("fgamma_i = " + fgamma_i + " fdelta = " + fdelta + " alpha = " + alpha);
                    if (rng.nextDouble() <= alpha) {
                        String s = delta.toString();
                        if (retval.indexOf(s) == -1) {
                            // Success: accept the individual and
                            // quit this k-loop
                            /// System.out.println("success, accepting.");
                            retval.add(s);
                            gamma_i = delta;
                            fgamma_i = fdelta;
                            break;
                        }
                        // System.out.println("fitness ok, alpha ok, but ind already in sample");
                    } else if (k > ntries) {
                        // Failure: ran out of tries. Accept the
                        // last individual only if it's not
                        // already in the sample
                        // System.out.println("failure, ran out of tries to find good ind");
                        String s = delta.toString();
                        if (retval.indexOf(s) == -1) {
                            // System.out.println("but accepting, ind wasn't in sample");
                            retval.add(s);
                        }
                        // Whether we accepted or not, quit this
                        // k-loop.
                        break;
                    } else {
                        // System.out.println("didn't accept because of alpha, looping.");
                    }
                }
            }
        }

        return retval;
    }

    // Get every possible node (up to a given depth). See the python
    // code python/RandomWalks/generate_trees.py for a more efficient
    // way to just print them.
    public ArrayList<String> sampleComplete() {

        ArrayList<String> retval = new ArrayList<String>();

        AllTrees at = new AllTrees(language);
        ArrayList<Node> allTrees = at.generateEntireSpace(maxDepth);
        for (Node t: allTrees) {
            String ts = t.toStringAsTree();
            retval.add(ts);
        }

        return retval;
    }


    // A central method: given two trees (as strings), calculate
    // syntactic and operator-based distances, and return them in a
    // hash.
    public HashMap<String, Double> getDistances(String s, String t) {

        HashMap<String, Double> retval = new HashMap<String, Double>();
        
        Tree sTree = new Tree(s);
        Tree tTree = new Tree(t);

		double logCutoff = 10e-44;

        // One-shot transition probability.
        double transProb = mutator.transitionProbability(sTree, tTree);
        retval.put("TP", transProb);

        // NCD
        double ncd = TreeDistance.getUniversalDistance(s.toString(), t.toString());
        retval.put("NCD", ncd);

        // FVD
        FeatureVectorDistance fvdObj = new
            FeatureVectorDistance(language.maxDepth, language.maxArity);
        double fvd = fvdObj.fvd(sTree, tTree);
        retval.put("FVD", fvd);
        FeatureVector sfv = new FeatureVector(sTree);
        FeatureVector tfv = new FeatureVector(tTree);
        retval.put("NodeCount", fvdObj.getNodeCount(sfv, tfv));
        retval.put("MinDepth", fvdObj.getMinDepth(sfv, tfv));
        retval.put("MeanDepth", fvdObj.getMeanDepth(sfv, tfv));
        retval.put("MaxDepth", fvdObj.getMaxDepth(sfv, tfv));
        retval.put("Symmetry", fvdObj.getSymmetry(sfv, tfv));
        retval.put("MeanFanout", fvdObj.getMeanFanout(sfv, tfv));
        retval.put("DiscreteMetric", fvdObj.getDiscreteMetricOnString(sfv, tfv));

        // TED FIXME could use the Augsten & Pawlik code instead?
        double ted = TreeDistance.getTreeEditDistance(sTree, tTree);
        retval.put("TED", ted);

        // TAD
        double tad[] = new double[6];
        for (int i = 0; i < 6; i++) {
            tad[i] = TreeDistance.getTreeAlignmentDistance(sTree, tTree, i);
            // System.out.println("TAD" + i + ": " + tad[i]);
            retval.put("TAD" + i, tad[i]);
        }

        // OVD
        double ovd = TreeDistance.getOverlapDistance(sTree, tTree);
        retval.put("OVD", ovd);

        return retval;
    }
        
    public static void main(String args[]) {

        if (args.length == 2 && args[0].equals("completeMatrices")) {
            int maxDepth = new Integer(args[1]);
            // write out the matrices of distances for the entire
            // space of given depth
            Sample sample = new Sample(maxDepth);
            sample.writeMatrices(sample.sampleComplete(),
                                 "depth_" + maxDepth,
                                 false);
            
        } else if (args.length == 2 && args[0].equals("uniformSampleMatrices")) {
            int maxDepth = new Integer(args[1]);
            // sample some individuals randomly and get distances
            // between them.
            Sample sample = new Sample(maxDepth);
            ArrayList<String> ofInterest = sample.sampleByGrow(10);
            sample.writeMatrices(ofInterest,
                                 "uniform_depth_" + maxDepth,
                                 true);
            
        } else if (args.length == 2 && args[0].equals("rwSampleMatrices")) {
            int maxDepth = new Integer(args[1]);
            // write out the matrices of distances for a sample from
            // the space of given depth, sampled by random walk.
            Sample sample = new Sample(maxDepth);
            sample.writeMatrices(sample.sampleRandomWalk(10, 100, 10),
                                 "../results/depth_" + maxDepth + "/rwSample",
                                 true);
            
        } else if (args.length == 2 && args[0].equals("mhSampleMatrices")) {
            int maxDepth = new Integer(args[1]);
            // write out the matrices of distances for a sample from
            // the space of given depth, sampled by
            // Metropolis-Hastings
            Sample sample = new Sample(maxDepth);
            sample.writeMatrices(sample.sampleMetropolisHastings(10, 20, 10),
                                 "mh_sample_depth_" + maxDepth,
                                 true);
            
        } else if (args.length == 2 && args[0].equals("randomWalking")) {
            int maxDepth = new Integer(args[1]);
            // sample some individuals randomly and estimate random
            // walk lengths between them by simulation.
            Sample sample = new Sample(maxDepth);
            ArrayList<String> ofInterest = sample.sampleByGrow(100);
            //sample.randomWalking(1000000000, ofInterest, 100);
            sample.randomWalking((int) (1298*1000), ofInterest, 100);

        } else if (args.length == 2 && args[0].equals("sampleForSuperNode")) {
            int maxDepth = new Integer(args[1]);
            // sample two individuals randomly, then get some near
            // neighbours of both, then calculate all pairwise
            // distances, then model the remainder of the graph as a
            // single supernode and calculate its distances as well.
            Sample sample = new Sample(maxDepth);
            String basename = ("../results/depth_" + maxDepth
                               + "/TP_supernode_estimates_RW/");
            sample.estimateTPWithSuperNode(basename, 50, 20);

        } else if (args.length == 2 && args[0].equals("sampleOneStepProbabilities")) {
            int maxDepth = new Integer(args[1]);
            // perform many mutations to estimate transition probabilities
            Sample sample = new Sample(maxDepth);
            sample.sampleOneStepProbabilities(sample.sampleComplete(),
                                              100,
                                              "../results/depth_" + maxDepth + "/TP_sampled.dat"
                                              );

        } else if (args.length == 2 && args[0].equals("testRandomWalk")) {
            int maxDepth = new Integer(args[1]);
            Sample sample = new Sample(maxDepth);
            for (String s: sample.sampleRandomWalk(10, 100, 10)) {
                System.out.println(s);
            }
            
        } else {
            System.out.println("Please read the source to see usage.");
        }
    }
}


