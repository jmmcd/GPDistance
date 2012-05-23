/*
 * TreeDistance.java
 *
 * Created on 3 November 2010
 *
 * Calculate feature-vector, tree-edit, tree-alignment, and normalised
 * compression distances on individuals' genotypes.
 *
 * @author James McDermott
 */

package TreeDistance;

import TreeDistance.TreeEditDistance.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

public class TreeDistance {

    // Normalised Compression Distance
    //////////////////////////////////

    // /**
    //  * @param i1, i2: GPIndividuals
    //  */
    // public static double getUniversalDistance(GPIndividual i1, GPIndividual i2) {
    //     return NormalisedCompressionDistance.ncd(i1.trees[0].toString(),
    //                                              i2.trees[0].toString());
    // }

    /**
     * @param s1, s2: Strings
     */
    public static double getUniversalDistance(String s1, String s2) {
        return NormalisedCompressionDistance.ncd(s1, s2);
    }

    /**
     * @param t1, t2: Trees
     */
    public static double getUniversalDistance(Tree t1, Tree t2) {
        return NormalisedCompressionDistance.ncd(t1, t2);
    }


    /**
     * @param t1, t2: Trees
     */
    public static double getFeatureVectorDistance(Tree t1, Tree t2) {
        FeatureVectorDistance fvd = new FeatureVectorDistance(2, 2);
        return fvd.fvd(t1, t2);
    }

    // Tree Alignment Distance
    //////////////////////////

    public static double getTreeAlignmentDistance(Tree i1, Tree i2) {
        // Discrete metric, depth weighting = 0.5
        TreeAlignmentDistance tads = new TreeAlignmentDistance(0.5);
        return tads.getDistance(i1, i2);
    }


    public static double getTreeAlignmentDistance(Tree i1, Tree i2, int variant) {
        TreeAlignmentDistance tads = null;
        switch (variant) {
        case 0:
            // Discrete metric, no depth weighting
            tads = new TreeAlignmentDistance(1.0);
            break;
        case 1:
            // Discrete metric, depth weighting = 0.5
            tads = new TreeAlignmentDistance(0.5);
            break;
        case 2:
            // Discrete metric, depth weighting = 0.1
            tads = new TreeAlignmentDistance(0.1);
            break;
        case 3:
            // Codes, no depth weighting
            tads = new TreeAlignmentDistance(makeSymbolicRegressionPrimitives(), 1.0);
            break;
        case 4:
            // Codes, depth weighting = 0.5
            tads = new TreeAlignmentDistance(makeSymbolicRegressionPrimitives(), 0.5);
            break;
        case 5:
            // Codes, depth weighting = 0.1
            tads = new TreeAlignmentDistance(makeSymbolicRegressionPrimitives(), 0.1);
            break;
        default:
            // Discrete metric, depth weighting = 0.5
            tads = new TreeAlignmentDistance(0.5);
        }
        return tads.getDistance(i1, i2);
    }

    public static List<String> makeSymbolicRegressionPrimitives() {
        String prims[] = {"+", "-", "*", "/", "x", "y"};
        return Arrays.asList(prims);
    }

    // /**
    //  * @param s1, s2: Strings to be interpreted as lisp-style
    //  * s-expressions representing trees.
    //  */
    // public static double getTreeAlignmentDistance(String s1, String s2) {
    //     GPTree gpt1 = GPTree.makeTreeFromSExpression(s1);
    //     GPTree gpt2 = GPTree.makeTreeFromSExpression(s2);
    //     TreeAlignmentDistance tads = new TreeAlignmentDistance();
    //     return tads.getDistance(gpt1, gpt2);
    // }



    // // Tree Edit Distance
    // /////////////////////

    // /**
    //  * @param s1, s2: Strings to be interpreted as lisp-style
    //  * s-expressions representing trees.
    //  */
    // public static double getTreeEditDistance(String s1, String s2) {
    //     CreateTreeHelper cth = new CreateTreeHelper();

    //     GPTree gpt1 = GPTree.makeTreeFromSExpression(s1);
    //     BasicTree bt1 = cth.makeTreeFromECJGPTree(gpt1);
    //     GPTree gpt2 = GPTree.makeTreeFromSExpression(s2);
    //     BasicTree bt2 = cth.makeTreeFromECJGPTree(gpt2);

    //     ComparisonZhangShasha treeCorrector = new ComparisonZhangShasha();
    //     OpsZhangShasha costs = new OpsZhangShasha();
    //     Transformation transform = treeCorrector.findDistance(bt1, bt2, costs);
    //     return transform.getCost();
    // }

    // /**
    //  * @param i1, i2: GPIndividuals
    //  */
    // public static double getTreeEditDistance(GPIndividual i1, GPIndividual i2) {
    //     CreateTreeHelper cth = new CreateTreeHelper();

    //     BasicTree bt1 = cth.makeTreeFromECJGPTree(i1.trees[0]);
    //     BasicTree bt2 = cth.makeTreeFromECJGPTree(i2.trees[1]);

    //     ComparisonZhangShasha treeCorrector = new ComparisonZhangShasha();
    //     OpsZhangShasha costs = new OpsZhangShasha();
    //     Transformation transform = treeCorrector.findDistance(bt1, bt2, costs);
    //     return transform.getCost();
    // }

    /**
     * @param i1, i2: Trees
     */
    public static double getTreeEditDistance(Tree i1, Tree i2) {
        CreateTreeHelper cth = new CreateTreeHelper();
        BasicTree bt1 = cth.makeBasicTreeFromTree(i1);
        BasicTree bt2 = cth.makeBasicTreeFromTree(i2);
        ComparisonZhangShasha treeCorrector = new ComparisonZhangShasha();
        OpsZhangShasha costs = new OpsZhangShasha();
        Transformation transform = treeCorrector.findDistance(bt1, bt2, costs);
        return transform.getCost();
    }


    /**
     * @param i1, i2: Trees
     */
    public static double getOverlapDistance(Tree i1, Tree i2) {
        return OverlapDistance.overlap(i1, i2);
    }


	public static void main(String args[]) {
        Tree t = new Tree("(* (* x x) y (/ (+ x))");
        Tree s = new Tree("(* (+ x y))");
        System.out.println("NCD: " + getUniversalDistance(t, s));
        System.out.println("FVD: " + getFeatureVectorDistance(t, s));
        System.out.println("TED: " + getTreeEditDistance(t, s));
        for (int i = 0; i < 6; i++) {
            System.out.println("TAD" + i + ": " + getTreeAlignmentDistance(t, s, i));
        }
        System.out.println("OVD: " + getOverlapDistance(t, s));
	}

}

