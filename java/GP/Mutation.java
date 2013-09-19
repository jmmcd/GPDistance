package GP;

import TreeDistance.*;

import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.ext.*;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.alg.KShortestPaths;

import java.util.*;
import java.io.PrintWriter;

import static java.lang.Math.*;

public class Mutation {


    Language language;
    int maxDepth;
    Random rng;
    int curID;

    /**
     * Construct the Mutation operator.
     *
     * @param _maxDepth maximum tree depth.
     *
     */
    public Mutation(Language _language, int _maxDepth, Random _rng) {
        language = _language;
        maxDepth = _maxDepth;
        rng = _rng;
        curID = 0;
    }

    /**
     * @param t The original Tree. It will not be altered. A mutated
     * version will be returned.
     */
    public Tree mutate(Tree t) {
		return mutate(t, true);
	}

    /**
     * @param t The original Tree. It will not be altered. A mutated
     * version will be returned.
     * @param allowRoot. Is mutation at the root allowed?
     */
    public Tree mutate(Tree t, boolean allowRoot) {
        Tree copy = t.clone();
        ArrayList<Node> nodes = copy.getRoot().depthFirstTraversal();
        int nnodes = nodes.size();
		int whichNode;
		if (allowRoot) {
			whichNode = rng.nextInt(nnodes);
		} else {
			whichNode = 1 + rng.nextInt(nnodes - 1);
		}
        Node n = nodes.get(whichNode);
        // System.out.println("Node picked is: " + whichNode + "; " + n.toString());
        int curDepth = n.getDepth();
        grow(n, maxDepth - curDepth);
        // System.out.println("New subtree there is: " + n.toStringAsTree());

        // relabel all nodes.
        int i = 0;
        for (Node node: copy.getRoot().depthFirstTraversal()) {
            node.id = i++;
        }
        return copy;
    }

    // Given a node, apply the grow algorithm to make it into a tree.
    public void grow(Node node, int depth) {
        // See node as a placeholder: we overwrite its label and
        // children, but we don't change its parent or id.
        node.children.clear();
        if (depth <= 0 || rng.nextDouble() < language.PT) {
            // make it a terminal
            node.label = language.chooseTerminal(rng);
        } else {
            // make it a function
            node.label = language.chooseFunction(rng);
            for (int i = 0; i < language.arities.get(node.label); i++) {
                Node newNode = new Node(node, curID++, "");
                node.children.add(newNode);
                grow(newNode, depth - 1);
            }
        }
    }


    // Assumes that t and s are proper trees, where IDs reflect
    // depth-first traversal and depths are calculated correctly.

    // This is wrong:
    // From (+ y (* x x))
    // To y
    // TP 0.0

    // From x
    // To y
    // TP 0.0

    // From (+ x y)
    // To y
    // TP 0.0    
    public float transitionProbability(Tree t, Tree s) {
        int F = language.F;
        int T = language.T;

        ArrayList<Integer> cutPointIDs = new ArrayList<Integer>();
        ArrayList<Node> requiredSubtrees = new ArrayList<Node>();
        ArrayList<Integer> cutPointDepths = new ArrayList<Integer>();

        ArrayList<Node> tNodes = t.getRoot().depthFirstTraversal();
        ArrayList<Node> sNodes = s.getRoot().depthFirstTraversal();

        ArrayList<Node> diffPts = differencePoints(t.getRoot(), s.getRoot());
        if (diffPts.size() == 0) {
            // if there are no differences, trees are identical and can cut
            // anywhere
            for (Node n: tNodes) {
                cutPointIDs.add(n.id);
                cutPointDepths.add(n.getDepth());
                requiredSubtrees.add(n);
            }
        } else {
            ArrayList<ArrayList<Integer>> ancestorPathsOfDiffPts
                = new ArrayList<ArrayList<Integer>>();
            // for each difference point in t, find a list of
            // ancestors by traversing up via parent links.
            for (Node diffPt: diffPts) {
                // System.out.println("diffPt: " + diffPt.toStringAsTree());
                ArrayList<Node> ancestorPath = diffPt.getAncestors();
                ArrayList<Integer> ancestorPathIDs = new ArrayList<Integer>();
                for (Node ancestor: ancestorPath) {
                    // System.out.println("ancestor of diffPt: " + ancestor);
                    // System.out.println("ancestorID of diffPt: " + ancestor.id);
                    ancestorPathIDs.add(ancestor.id);
                }
                ancestorPathsOfDiffPts.add(ancestorPathIDs);
            }
            // can cut at any of the common ancestors.
            for (Integer n: findCommon(ancestorPathsOfDiffPts)) {
                // System.out.println("Common ancestor: " + n);
                cutPointIDs.add(n);
            }
            // key point: because the possible cutpoints will be a
            // single path from root downward, under which *all
            // differences* are contained, a DF traversal of t and s
            // will give the same ids for them. so iterate through s
            // and when we find a cutpoint id, take the subtree from
            // there as the subtree that would need to be generated to
            // make s, if t was cut there. Once we have a
            // requiredSubtree for each ID, can quit the loop.
            for (int i = 0;
                 i < sNodes.size()
                     && requiredSubtrees.size() < cutPointIDs.size();
                 i++) {

                Node n = sNodes.get(i);
                if (cutPointIDs.contains(n.id)) {
                    cutPointDepths.add(n.getDepth());
                    requiredSubtrees.add(n);
                }
            }
        }

        // Then just need to sum over 1/nnodes * prob(requiredSubtree).
        float retval = 0.0f;
        int nnodes = tNodes.size();
        // System.out.println("nnodes = " + nnodes);
        for (int i = 0; i < requiredSubtrees.size(); i++) {

            // debug stuff
            // System.out.println("cutPoint: " + cutPointIDs.get(i));
            // System.out.println("cutPointDepth: " + cutPointDepths.get(i));
            // System.out.println("requiredSubtree:"
            //                    + requiredSubtrees.get(i).toStringAsTree());

            float sgp = subtreeGrowProb(requiredSubtrees.get(i),
                                        maxDepth - cutPointDepths.get(i));
            if (sgp < 0.0f) {
                // System.out.println("trans prob retval = " + 0.0f);
                return 0.0f;
            }
            retval += (1.0f / nnodes) * sgp;
        }
        // System.out.println("trans prob retval = " + retval);
        return retval;
    }



    // "best" cut point is the lowest common ancestor, ie that with
    // greatest depth. But we're iterating over all possible cutpoints
    // anyway.

    // Take the union of the ancestor-lists
    public static ArrayList<Integer>
        findCommon(ArrayList<ArrayList<Integer>> inputs) {

        ArrayList<Integer> current = inputs.get(0);

        for (int j = 1; j < inputs.size(); j++) {
            current.retainAll(inputs.get(j));
        }
        return current;
    }

    // perform depth-first traversal in lock-step at each difference,
    // add the node to a list of differences and return from that
    // subtree. problem: will the nodes/ids be identical between the
    // two trees? no: save them as identified in the first tree.
    // remember, return as soon as we see a difference.
    public static ArrayList<Node> differencePoints(Node t, Node s) {
        ArrayList<Node> retval = new ArrayList<Node>();
        if (!t.label.equals(s.label)) {
            retval.add(t);
        } else {
            if (t.children.size() > 0) {
                for (int i = 0; i < t.children.size(); i++) {
                    retval.addAll(differencePoints(t.children.get(i),
                                                   s.children.get(i)));
                }
            }
        }
        return retval;
    }


    // Probability of generating a particular subtree using
    // the grow method.
    public float subtreeGrowProb(Node n, int maxDepth) {

        int T = language.T;
        int F = language.F;
        // PT is the probability of choosing a terminal
        float PT = language.PT;

        if (maxDepth <= 0 && n.children.size() > 0) {
            // Something's gone wrong: this tree can't be created
            return -1.0f;
        }

        if (n.children.size() == 0) {
            // The subtree is a single terminal

            if (maxDepth == 0) {
                // we MUST choose a terminal
                // System.out.println("subtreeGrowProb for subtree "
                //                    + n.toStringAsTree() + " " + 1.0f / T);
                return 1.0f / T;
            } else {

                // we don't HAVE to choose a terminal, but we do.
            
                // (1 / T) is the probability of choosing correct
                // terminal, given that we choose a terminal
                // System.out.println("subtreeGrowProb for subtree "
                //                    + n.toStringAsTree() + " " + PT / T);
                return PT / T;
            }
        }

        // The subtree's root is a function.  (1 - PT) is the
        // probability of choosing a function. (1 / F) is the
        // probability of choosing correct function, given that we
        // choose a function.
        float retval = (1.0f - PT) / F;
        for (Node child: n.children) {
            // recurse for each child.
            float sgp = subtreeGrowProb(child, maxDepth - 1);
            if (sgp < 0.0f) {
                return -1.0f;
            }
            retval *= sgp;
        }
        // System.out.println("subtreeGrowProb for subtree "
        //                    + n.toStringAsTree() + " " + retval);
        return retval;
    }

    // For full trees, every tree is equally likely. So ignore the
    // passed-in node.
    public float subtreeFullProb(Node n, int maxDepth) {
        return 1.0f / language.countFullTrees();
    }


    // Mutate an individual lots of times and keep the good ones.
    // Unused for now.
	public ArrayList<Tree> mutateAndSelect(Tree ind,
                                           int nMutations,
                                           int selectionSize) {

        // Generate individuals by mutation from ind.
        ArrayList<Tree> pop = new ArrayList<Tree>();
		for (int i = 0; i < nMutations; i++) {
			pop.add(mutate(ind));
		}

        if (nMutations == selectionSize) {
            // No need to sort or select.
            return pop;
        }

		// Sort in ascending order, that is best-first.
		Collections.sort(pop, new individualFitnessComparator());

		// Discard the bad ones.
		pop.subList(selectionSize, pop.size()).clear();
		return pop;
	}

    // Mutate each individual in a population lots of times and keep
    // the good ones. Unused for now.
	public ArrayList<Tree> mutatePopAndSelect(ArrayList<Tree> pop,
                                              int nMutations,
                                              int selectionSize) {

        // Generate individuals by mutation from pop.
        ArrayList<Tree> newpop = new ArrayList<Tree>();
        for (Tree t: pop) {
            for (int i = 0; i < nMutations; i++) {
                newpop.add(mutate(pop.get(i)));
            }
        }

		// Sort in ascending order, that is best-first.
		Collections.sort(newpop, new individualFitnessComparator());

		// Discard the bad ones.
		newpop.subList(selectionSize, newpop.size()).clear();
		return newpop;
	}

    // How to compare individuals' fitness for sorting methods. FIXME
    // there should be a boolean somewhere for minimise/maximise
    // fitness. For now we assume fitness is always minimised.
    public class individualFitnessComparator implements Comparator {

        public individualFitnessComparator() {
        }

        public int compare(Object o1, Object o2) {
            float f1 = (float) ((Tree) o1).fitness();
            float f2 = (float) ((Tree) o2).fitness();
            if (f1 < f2) {
                return -1;
            } else if (f1 > f2) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    // Test: make two trees and calculate the transition probability.
    public static void main(String args[]) {
        if (args.length != 2) {
            System.out.println("Usage: <tree1> <tree2>");
            System.exit(1);
        }
        Tree t = new Tree(args[0]);
        Tree s = new Tree(args[1]);
        int maxDepth = 2;
        Mutation mutator = new Mutation(new Language(maxDepth),
                                        maxDepth, new Random());
        System.out.println(mutator.transitionProbability(t, s));
    }
}


