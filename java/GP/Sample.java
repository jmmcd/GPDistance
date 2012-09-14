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

    // Estimate the length of a random walk by simulation. This
    // function doesn't collect a sample of individuals. It just
    // performs walks and saves the number of steps between pairs of
    // individuals. So it aims to estimate FMPT by simulation rather
    // than exact methods on a sample from the space.
    public HashMap<String, HashMap<String, Double>>
        randomWalking(ArrayList<String> ofInterest) {

        // lastOccur holds the individuals of interest as strings and
        // the number of steps since their last occurence, as
        // integers. These are created first.  Then when we're
        // walking, anytime we encounter one of these we learn
        // something about it.
        HashMap<String, Integer> lastOccur =
            new HashMap<String, Integer>();

        // walkLengths holds pairs of individuals of interest
        // (strings) and a list of samples of walkLengths between
        // them.
        HashMap<String, HashMap<String, ArrayList<Integer>>> walkLengths =
            new HashMap<String, HashMap<String, ArrayList<Integer>>>();

        // set up lastOccur
        for (String s: ofInterest) {
            lastOccur.put(s, -1);
        }

        // set up walkLengths
        for (String s: lastOccur.keySet()) {
            HashMap<String, ArrayList<Integer>> tmp =
                new HashMap<String, ArrayList<Integer>>();
            for (String t: lastOccur.keySet()) {
                tmp.put(t, new ArrayList<Integer>());
            }
            walkLengths.put(s, tmp);
        }

        // start with a random node
        Tree current = new Tree("x");
        mutator.grow(current.getRoot(), maxDepth);
        String cs = current.toString();
        
        // perform a random walk
        for (int i = 0; i < 1000; i++) {

            Integer csLastOccur = lastOccur.get(cs);
            if (csLastOccur == null) {

                // It's not one of our individuals of interest
                // +1 to all non-negative entries in lastOccur
                for (String t: lastOccur.keySet()) {
                    int tLastOccur = lastOccur.get(t);
                    if (tLastOccur != -1) {
                        lastOccur.put(t, lastOccur.get(t) + 1);
                    }
                }
                
            } else {
                
                // It's one of our individuals of interest: for all
                // non-negative t in lastOccur, save a sample (t, s)
                // in walkLengths
                for (String t: lastOccur.keySet()) {
                    int tLastOccur = lastOccur.get(t);
                    if (tLastOccur != -1) {
                        walkLengths.get(t).get(cs).add(tLastOccur);
                    }
                }
                // Mark the last occurrence as now.
                lastOccur.put(cs, 0);
            }

            current = mutator.mutate(current);
            cs = current.toString();
        }

        // Get the mean walk lengths and return them. FIXME should
        // possibly consider assuming a distribution, fitting that
        // distribution, and returning the mean?
        HashMap<String, HashMap<String, Double>> meanWalkLengths =
            new HashMap<String, HashMap<String, Double>>();
        for (String s: lastOccur.keySet()) {
            HashMap<String, Double> tmp =
                new HashMap<String, Double>();
            for (String t: lastOccur.keySet()) {
                tmp.put(t, mean(walkLengths.get(t).get(s)));
            }
            meanWalkLengths.put(s, tmp);
        }
        return meanWalkLengths;
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


    // A central method: give two trees (as strings), calculate
    // syntactic and operator-based distances, and return them
    // in a hash.
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
            
        } else if (args.length == 2 && args[0].equals("randomWalkingMatrices")) {
            int maxDepth = new Integer(args[1]);
            // sample some individuals randomly and estimate random
            // walk lengths between them by simulation.
            Sample sample = new Sample(maxDepth);
            ArrayList<String> ofInterest = sample.sampleUniform(10);
            HashMap<String, HashMap<String, Double>>
                hshsali = sample.randomWalking(ofInterest);
            sample.writeMatrices(ofInterest,
                                 "random_walking_" + maxDepth,
                                 true, hshsali);
            
        } else {
            System.out.println("Please read the source to see usage.");
        }
    }
}
