package GP;

import TreeDistance.*;

import static java.lang.Math.*;
import java.util.*;

/*
 * AllTrees: Make a list of all trees up to a given depth.
 *
 * We use a bottom-up approach. Build small trees and recombine them.
 * For this we need external storage -- the treesOfDepth and
 * treesOfDepthLE variables. This avoids re-calculating small trees
 * zillions of times over. The Node "x", for example, is reused a lot
 * of times, rather than having a new Node every time it's required.
 * One disadvantage of this approach is that each tree is not truly
 * instantiated separately, so nodes don't know their real parents,
 * for example. The user of the class can clone any tree in order to
 * get a proper tree with correct ids, parent links and so on.
 *
 */

public class AllTrees {

    // Trees of a given depth.
    public ArrayList<ArrayList<Node>> treesOfDepth;
    // Trees of a given depth or less.
    public ArrayList<ArrayList<Node>> treesOfDepthLE;
    public Node dummy;
    public Language language;

    public AllTrees(Language _language) {

        treesOfDepth = new ArrayList<ArrayList<Node>>();
        treesOfDepthLE = new ArrayList<ArrayList<Node>>();

        language = _language;
    }

    // Top-level method.
    public ArrayList<Node> generateEntireSpace(int maxDepth) {

        for (int i = 0; i <= maxDepth; i++) {
            treesOfDepth.add(generateTreesOfDepth(i));
            treesOfDepthLE.add(new ArrayList<Node>());
            for (int j = 0; j <= i; j++) {
                treesOfDepthLE.get(i).addAll(treesOfDepth.get(j));
            }
        }
        return treesOfDepthLE.get(maxDepth);
    }

    // Calculate trees of a given depth
    public ArrayList<Node> generateTreesOfDepth(int depth) {

        ArrayList<Node> retval = new ArrayList<Node>();

        if (depth == 0) {
            // It's just leaves
            for (String t: language.terminals) {
                Node n = new Node(dummy, 0, t);
                retval.add(n);
            }
            return retval;
        }

        // Every possible function combined with every possible
        // assignment of smaller subtrees to function's children.
        for (String f: language.functions) {
            for (ArrayList<Node> a:
                     generateAssignment(treesOfDepthLE.get(depth - 1),
                                        language.arities.get(f))) {
                Node n = new Node(dummy, 0, f);
                n.children = a;
                if (n.cloneAsTree().getMaxDepth() == depth) {
                    // Add n only if at least one element of a was
                    // large enough to make n have the right depth.
                    retval.add(n);
                }
            }
        }
        return retval;
    }

    // Given a list of things from which we must pick n times with
    // replacement (such that order is important), generate all the
    // possible outcome sequences. Eg vals = [a, b, c], n = 2 returns:
    // [aa, ab, ac, ba, bb, bc, ca, cb, cc].
    public ArrayList<ArrayList<Node>>
        generateAssignment(ArrayList<Node> vals, int n) {

        ArrayList<ArrayList<Node>> retval = new ArrayList<ArrayList<Node>>();
        if (n == 1) {
            for (Node val: vals) {
                ArrayList<Node> tmp = new ArrayList<Node>();
                tmp.add(val);
                retval.add(tmp);
            }
            return retval;
        } else {
            for (Node val: vals) {
                for (ArrayList<Node> rest: generateAssignment(vals, n-1)) {
                    ArrayList<Node> tmp = new ArrayList<Node>();
                    tmp.add(val);
                    tmp.addAll(rest);
                    retval.add(tmp);
                }
            }
            return retval;
        }
    }

    public static void main(String args[]) {
        int maxDepth = Integer.parseInt(args[0]);
        AllTrees generator = new AllTrees(new Language(maxDepth));
        ArrayList<Node> space = generator.generateEntireSpace(maxDepth);
        HashMap<String, Boolean> names = new HashMap<String,Boolean>();
        for (Node n: space) {
            String name = n.toStringAsTree();
            System.out.println(name);
            if (names.get(name) != null) {
                System.out.println("found a duplicate! " + name);
            }
            names.put(n.toStringAsTree(), true);
        }
        System.out.println("Number of trees of depth <= " + 
                           maxDepth + ": " + space.size());
    }
}
