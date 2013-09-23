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

    // Given a Hashmap of Hashmaps of values, write out a matrix.
    public void writeMatrix(ArrayList<String> trees,
                            HashMap<String, HashMap<String, Double>> vals,
                            String filename) {
        try {
            // Open file
            FileWriter fw = new FileWriter(filename);

            // for every source tree
            for (String s: trees) {
                HashMap<String, Double> vals_s = vals.get(s);
                // for every destination tree
                for (String t: trees) {
                    fw.write(vals_s.get(t) + " ");
                }
                fw.write("\n");
            }
            
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
            
        
    

    // Given a list of trees, write out a matrix of distances between
    // all pairs, for many different types of distance. Can pass in
    // results from a random-walking simulation as mfpte (pass null
    // otherwise).
    public void writeMatrices(ArrayList<String> trees,
                              String codename,
                              boolean writeNonNormalised,
                              HashMap<String, HashMap<String, Double>> mfpte
                              ) {

        ArrayList<String> distanceNames = new ArrayList<String>();
        String _distanceNames[] = {
            "TP", "NCD", "FVD",
            "NodeCount", "MinDepth", "MeanDepth", "MaxDepth",
            "Symmetry", "MeanFanout", "DiscreteMetric",
            "TED",
            "TAD0", "TAD1", "TAD2", "TAD3", "TAD4", "TAD5",
            "OVD"
        };

        for (String d: _distanceNames) {
            distanceNames.add(d);
        }

        if (mfpte != null) {
            distanceNames.add("MFPTE");
        }        

        HashMap<String, FileWriter> files = new HashMap<String, FileWriter>();
        try {
            // Open files
            for (String distance: distanceNames) {

                // write to ../results/<codename>/<distance>.dat
                String dirname = "../results/" + codename + "/";
                (new File(dirname)).mkdirs();
                String filename;
                if (distance.equals("TP") && writeNonNormalised) {
                    // Since we're sampling from the space, the
                    // transition probabilities won't sum to 1. So
                    // save to _nonnormalised.
                    filename = dirname + "TP_nonnormalised.dat";
                } else {
                    filename = dirname + distance + ".dat";
                }
                files.put(distance, new FileWriter(filename)); 
            }

            // for every source tree
            for (String s: trees) {
                // for every destination tree
                for (String t: trees) {

                    // Get distances...
                    HashMap<String, Double> distances = getDistances(s, t);

                    // ...hack in the mfpte value if there is one,
                    // else zero...
                    if (mfpte != null) {
                        distances.put("MFPTE", mfpte.get(s).get(t));
                    }

                    // ... and write them out
                    for (String distance: distanceNames) {
                        files.get(distance).write(distances.get(distance) + " ");
                    }
                }
                // write a newline
                for (String distance: distanceNames) {
                    files.get(distance).write("\n");
                }
            }

            // close files
            for (String distance: distanceNames) {
                files.get(distance).close();
            }
        
        } catch (IOException e) {
            e.printStackTrace();
        }
        
    }


    public ArrayList<String> sampleUniform(int n) {

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


    // TODO alternative methods of deriving an estimate of FMPT for
    // pairs of nodes to confirm the FMPT calculations?
    

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
                float[] counts = new float[trees.size()];
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
                float sum = 0.0f;
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
        
            float meanRowsum = 0.0f; // sum of inward TPs to supernode over all rows but last
            for (int i = 0; i < L; i++) {
                float rowsum = 0.0f;
                for (int j = 0; j < L; j++) {
                    float tp = mutator.transitionProbability(new Tree(sample.get(i)),
                                                             new Tree(sample.get(j)));
                    rowsum += tp;
                    tpFile.write(tp + " ");

                    // System.out.println("From " + sample.get(i));
                    // System.out.println("To " + sample.get(j));
                    // System.out.println("TP " + tp);

                }

                // write the inward TP for supernode
                tpFile.write(1.0f - rowsum + "");
                meanRowsum += (1.0f - rowsum);

                // System.out.println("rowsum " + rowsum);
                tpFile.write("\n");
            }

            meanRowsum /= L;

            // TP(x, S) = 1 - sum(TP(x, y)) for all y = 1 - rowsum
            
            // assume TP(S, S) = mean(rowsum)

            // TP(S, x) = (1 - mean(rowsum)) / L -- share out remaining probability

            for (int i = 0; i < L; i++) {
                tpFile.write((1.0f - meanRowsum) / L + " ");
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
    // than exact methods on a sample from the space. FIXME should aim
    // to test this on a known space, eg the 1298-space.

    // We only need to save the last occurrence of each of the
    // individuals of interest. Don't need the complete history, which
    // may be in the billions. Each time we hit an individual of
    // interest s, we can save a sampled walklength(t, s) =
    // lastOccur(s) - lastOccur(t) if we have previously seen t. The
    // picture is: others -> t -> others -> s. Note that it can happen
    // that we see t multiple times before seeing s: t0 -> others ->
    // t1 -> others -> t2 -> others -> s. In that case we should save
    // (s-t0), but we should not save (s-t1) or (s-t2). We can also
    // see s multiple times: t -> others -> s0 -> others -> s1 ->
    // others -> s2. In that case we save (s0-t), but we should not
    // save (s1-t) or (s2-t). to achieve this, the current
    // implementation below (lastOccur = new HashMap<String, Long>())
    // is insufficient -- it will not avoid saving the samples we
    // shouldn't save. can we use arrays to store? maybe put a limit
    // of say 100 samples and store in a [][][]?

    // we need to know size of space. the expected ELRW is size/2. if
    // that's much bigger than the largest feasible walk, then we
    // can't get good sampling, only biased sampling
    
    public void randomWalking(ArrayList<String> ofInterest,
                              String filename) {

        // lastOccur holds the individuals of interest as strings and
        // the number of steps since their last occurence, as
        // integers. These are created first.  Then when we're
        // walking, anytime we encounter one of these we learn
        // something about it.

        HashMap<String, Long> lastOccur =
            new HashMap<String, Long>();

        // walkLengths holds pairs of individuals of interest
        // (strings) and a list of samples of walkLengths between
        // them.
        HashMap<String, HashMap<String, ArrayList<Long>>> walkLengths =
            new HashMap<String, HashMap<String, ArrayList<Long>>>();

        // set up lastOccur
        for (String s: ofInterest) {
            lastOccur.put(s, -1L);
        }

        // set up walkLengths
        for (String s: ofInterest) {
            HashMap<String, ArrayList<Long>> tmp =
                new HashMap<String, ArrayList<Long>>();
            for (String t: ofInterest) {
                tmp.put(t, new ArrayList<Long>());
            }
            walkLengths.put(s, tmp);
        }

        // start at one of the individuals of interest
        String cs = ofInterest.get(0);
        Tree current = new Tree(cs);

        long lim = 10000;
        // perform a random walk
        for (long i = 0; i < lim; i += 1) {

            Long csLastOccur = lastOccur.get(cs);
            if (csLastOccur == null) {

                System.out.println(cs + " not of interest");
                // It's not one of our individuals of interest
                // +1 to all non-negative entries in lastOccur
                for (String t: lastOccur.keySet()) {
                    long tLastOccur = lastOccur.get(t);
                    if (tLastOccur != -1L) {
                        lastOccur.put(t, lastOccur.get(t) + 1);
                    }
                }
                for (String t: lastOccur.keySet()) {
                    System.out.print(lastOccur.get(t) + " ");
                }
                System.out.println("\n");

            } else {
                
                // It's one of our individuals of interest: for all
                // non-negative t in lastOccur, save a sample (t, s)
                // in walkLengths
                for (String t: lastOccur.keySet()) {
                    long tLastOccur = lastOccur.get(t);
                    if (tLastOccur != -1) {
                        // +1 because if we mutate immediately to
                        // self, walk length should be 1 (not 0).
                        walkLengths.get(t).get(cs).add(tLastOccur + 1);
                        System.out.println(cs + " of interest: saving " + (tLastOccur+1));

                    }
                }
                // Mark the last occurrence as now.
                lastOccur.put(cs, 0L);

                for (String t: lastOccur.keySet()) {
                    System.out.print(lastOccur.get(t) + " ");
                }
                System.out.println("\n");

            }

            current = mutator.mutate(current);
            cs = current.toString();

            if (i % 10000 == 0) {
                System.out.println("" + i + " of " + lim + " mutations done");
            }
        }


        try {
            // Open file
            FileWriter fw = new FileWriter(filename);
        
            // Write out the samples for each pair, one pair per line, in
            // the natural order
            for (String t: lastOccur.keySet()) {
                HashMap<String, ArrayList<Long>> tmp = walkLengths.get(t);
                for (String s: lastOccur.keySet()) {
                    fw.write(t + ": ");
                    fw.write(s + ": ");
                    for (Long i: tmp.get(s)) {
                        fw.write(i + " ");
                    }
                    fw.write("\n");
                }
            }
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
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
                                 false, null);
            
        } else if (args.length == 2 && args[0].equals("uniformSampleMatrices")) {
            int maxDepth = new Integer(args[1]);
            // sample some individuals randomly and get distances
            // between them.
            Sample sample = new Sample(maxDepth);
            ArrayList<String> ofInterest = sample.sampleUniform(10);
            sample.writeMatrices(ofInterest,
                                 "uniform_depth_" + maxDepth,
                                 true, null);
            
        } else if (args.length == 2 && args[0].equals("rwSampleMatrices")) {
            int maxDepth = new Integer(args[1]);
            // write out the matrices of distances for a sample from
            // the space of given depth, sampled by random walk.
            Sample sample = new Sample(maxDepth);
            sample.writeMatrices(sample.sampleRandomWalk(10, 20, 10),
                                 "rw_sample_depth_" + maxDepth,
                                 true, null);
            
        } else if (args.length == 2 && args[0].equals("mhSampleMatrices")) {
            int maxDepth = new Integer(args[1]);
            // write out the matrices of distances for a sample from
            // the space of given depth, sampled by
            // Metropolis-Hastings
            Sample sample = new Sample(maxDepth);
            sample.writeMatrices(sample.sampleMetropolisHastings(10, 20, 10),
                                 "mh_sample_depth_" + maxDepth,
                                 true, null);
            
        } else if (args.length == 2 && args[0].equals("randomWalking")) {
            int maxDepth = new Integer(args[1]);
            // sample some individuals randomly and estimate random
            // walk lengths between them by simulation.
            Sample sample = new Sample(maxDepth);
            ArrayList<String> ofInterest = sample.sampleUniform(4);
            sample.randomWalking(ofInterest, "../results/depth_" + maxDepth
                                 + "/MFPT_random_walking_samples.dat");

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


