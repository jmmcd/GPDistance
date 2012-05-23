/*
 * OverlapDistance.java
 *
 * Created on June 9 2011
 *
 */
package TreeDistance;

import java.util.*;
import jscheme.JScheme;
import static java.lang.Math.*;

/**
 * Calculate Lee Spector's overlap distance between two trees.
 *
 * From http://push.i3ci.hampshire.edu/2011/06/07/overlap/
 *
 *
 * Overlap is defined in terms of the collections of the items
 * contained in each of the arguments, including nested items at all
 * levels. For example, if x is (a (b c)) then Items(x) will contain
 * the following five items: a, b, c, (b c) and (a (b c)). Within this
 * context we define:
 *
 * Overlap(x, y) = the maximal number of simultaneous, disjoint
 * pairings by identity of items in Items(x) with items in Items(y),
 * divided by the size of the larger of Items(x) and Items(y).
 *
 *
 * Overlap as defined by Lee Spector ranges from 0 (entirely distinct)
 * to 1 (entirely similar). We invert this to make it look like a
 * distance, returning 1 - overlap, so a similar pair gives a small
 * number and a dissimilar pair a larger one.
 *
 * Java implementation only by jmmcd.
 *
 * @author James McDermott
 */
public class OverlapDistance {

    // Two entry points: overlap(Tree, Tree) and overlap(Node, Node).
    // Use whichever suits.

    // An entry point. Calculate overlap.
    public static double overlap(Tree x, Tree y) {
        return overlap(x.getRoot(), y.getRoot());
    }

    // An entry point. Calculate overlap.
    public static double overlap(Node x, Node y) {
        int overlapCount = 0;

        // Calculate all the items in X and Y
        ArrayList<String> itemsX = allItems(x);
        ArrayList<String> itemsY = allItems(y);

        // Calculate their frequencies
        HashMap<String, Integer> freqsX = frequencies(itemsX);
        HashMap<String, Integer> freqsY = frequencies(itemsY);

        // For any key which occurs in both frequency maps, add to the
        // overlap-count.
        for (String key: freqsX.keySet()) {
            if (freqsY.containsKey(key)) {
                // this key is in both, so contributes to overlapCount.
                overlapCount += min(freqsX.get(key), freqsY.get(key));
            }
        }

        // overlap is a ratio between overlapCount and size.
        double overlap = overlapCount / (double) max(itemsX.size(), itemsY.size());

        // we invert the sense to make it look like a distance.
        return 1.0 - overlap;
    }

    // Get a list of all the nodes and subtrees descended from a given
    // node.
    private static ArrayList<String> allItems(Node x) {
        ArrayList<String> retval = new ArrayList<String>();
        if (x.children.size() > 0) {
            retval.add(x.toStringAsTree());
        }
        retval.add(x.toString());
        for (Node n: x.children) {
            retval.addAll(allItems(n));
        }
        return retval;
    }

    // Get a map from items to a count of their occurences in the
    // input list.
    private static HashMap<String, Integer> frequencies(ArrayList<String> a) {
        HashMap<String, Integer> retval = new HashMap<String, Integer>();
        for (String item: a) {
            if (retval.containsKey(item)) {
                retval.put(item, retval.get(item) + 1);
            } else {
                retval.put(item, 1);
            }
        }
        return retval;
    }


    public static void main(String []args) {

        // These test cases come from the same web-page
        // mentioned above.
        String tests = "" +
            "a , () --- float overlap: 0.0\n" +
            "a , a --- float overlap: 1.0\n" +
            "a , b --- float overlap: 0.0\n" +
            "a , (a b) --- float overlap: 0.33333334\n" +
            "a , (a a) --- float overlap: 0.33333334\n" +
            "(a) , (a b) --- float overlap: 0.33333334\n" +
            "(a b) , (a c) --- float overlap: 0.33333334\n" +
            "(a b) , (a b c) --- float overlap: 0.5\n" +
            "(a b c) , (a b d) --- float overlap: 0.5\n" +
            "(a b c d) , (a b c e) --- float overlap: 0.6\n" +
            "(a b c d) , (d c b a) --- float overlap: 0.8\n" +
            "(a b c d e f g) , (a b c d e f h) --- float overlap: 0.75\n" +
            "(a b) , (a b c d e f) --- float overlap: 0.2857143\n" +
            "(a (b (c))) , (a (b (c))) --- float overlap: 1.0\n" +
            "(a (b (c))) , (a (b (x))) --- float overlap: 0.33333334\n" +
            "(a (b (c))) , (a (x (c))) --- float overlap: 0.5\n" +
            "(a (b (c))) , (x (b (c))) --- float overlap: 0.6666667\n" +
            "(a (b c) (d (e f))) , ((d (e f)) a) --- float overlap: 0.6\n" +
            "(a (b c) (d (e f))) , (a a (b c) (d (e f))) --- float overlap: 0.8181818";

        for (String test: tests.split("\n")) {
            String [] parts = test.split("---");
            String [] trees = parts[0].split(",");
            String [] floats = parts[1].split(":");
            float v = Float.parseFloat(floats[1]);
            System.out.println("Tree A: " + trees[0] + "; Tree B: " + trees[1] + "; overlap: " + v);
            // Test allItems()
            Tree x = new Tree(trees[0]);
            for (String item: allItems(x.getRoot())) {
                System.out.println(item);
            }

            System.out.println("My result: " + overlap(new Tree(trees[0]), new Tree(trees[1])));

        }
    }
}


