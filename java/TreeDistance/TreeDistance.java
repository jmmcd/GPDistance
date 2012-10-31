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

/*
 * RTED is COPYRIGHT 2012 Mateusz Pawlik and Nikolaus Augsten
 * 
 * RTED 1.1 downloaded from
 * [http://www.inf.unibz.it/dis/projects/tree-distance-repository]
 *
 * Licence (from RTED homepage): Our program is free to redistribute
 * and/or modify under the terms of the GNU Affero General Public
 * License [http://www.gnu.org/licenses/].
 *
 * Citing RTED (from RTED homepage): If you want to refer to RTED in a
 * publication, please cite the following PVLDB paper.
 * 
 * M. Pawlik and N. Augsten. RTED: A Robust Algorithm for the Tree
 * Edit Distance. PVLDB 5(4):334–345, 2011.
 */

// util and distance are the unfortunate package names used by the
// RTED library.
import util.LblTree;
import distance.RTED_InfoTree_Opt;


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


    /**
     * @param i1, i2: Trees
     */
    public static double getOverlapDistance(Tree i1, Tree i2) {
        return OverlapDistance.overlap(i1, i2);
    }


    /**
     * @param n: a Node. On first call, n is the root of a Tree. Then
     * we recurse with n equal to each of the children.
     *
     * @return a tree in the format used by the RTED tree distance
     * library.
     */ 
    public static LblTree lblTreeFromTree(Node n) {
        int treeID = 0;
        LblTree node = new LblTree(n.toString(), treeID);
		for (Node child: n.children) {
			node.add(lblTreeFromTree(child));
		}
		return node;
    }


    /**
     * @param t1, t2: Trees. We'll convert them to LblTree format as
     * used by RTED.
     *
     * This implementation uses the RTED library to calculate tree
     * edit distance. It should be preferred because the other
     * implementation could be flawed, see
     * [https://github.com/irskep/zhang-shasha]
     *
     * @return the tree-edit distance, a nonnegative double value.
     */
    public static double getTreeEditDistance(Tree t1, Tree t2) {
        LblTree lt1 = lblTreeFromTree(t1.getRoot());
        LblTree lt2 = lblTreeFromTree(t2.getRoot());
        // 1, 1, 1 indicate the insert, delete, edit costs
		RTED_InfoTree_Opt rted = new RTED_InfoTree_Opt(1, 1, 1);
        rted.init(lt1, lt2);
        rted.computeOptimalStrategy();
        return rted.nonNormalizedTreeDist();
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
