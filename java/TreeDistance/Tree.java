/*
 * Tree.java
 *
 * Created November 3 2010
 *
 */

package TreeDistance;

import jscheme.JScheme;

import java.util.ArrayList;
import java.util.StringTokenizer;
import static java.lang.Math.*;


public class Tree {

    // holder is a node whose only child is the root
    Node holder;
    int currID;
    Node currNode;
    int nodeCount;
    int maxDepth;
    int minDepth;
    float meanDepth;
    float symmetry;
    float meanFanout;
    String program = null;

    protected final JScheme js;

    // Constructor from a string.
    public Tree(String input) {
        nodeCount = -1;
        maxDepth = -1;
        minDepth = 1000;
        meanDepth = -1.0f;
        symmetry = -1.0f;
        meanFanout = -1.0f;

        js = new JScheme();
        String pdiv = "(define p (lambda (a b) (if (< (abs b) 0.00001) a (/ a b))))";
        js.eval(pdiv);

        // Make sure the string is tokenizable
        // FIXME allow other delimiters?
        input = input.replace("(", " ( ");
        input = input.replace("[", " [ ");
        input = input.replace(")", " ) ");
        input = input.replace("]", " ] ");

        currID = 0;
        holder = new Node(null, currID++, "holder");
        StringTokenizer st = new StringTokenizer(input);
        parseString(holder, st);
    }

	// FIXME make a random constructor
	// // Generate a random individual
	// public Tree(int maxDepth) {
    //     nodeCount = -1;
    //     maxDepth = -1;
    //     minDepth = 1000;
    //     meanDepth = -1.0f;
    //     symmetry = -1.0f;
    //     meanFanout = -1.0f;

	// 	js = new JScheme();
    //     String pdiv = "(define p (lambda (a b) (if (< (abs b) 0.00001) a (/ a b))))";
    //     js.eval(pdiv);

	// 	currID = 0;
    //     holder = new Node(null, currID++, "holder");


    public Tree clone() {
        // hehe, the old serialise-deserialise trick
        return new Tree(toString());
    }

    public Node getRoot() {
        return holder.children.get(0);
    }

    public void calculateFeatures() {
        int depthSum = 0;
        int nChildren = 0;
        int nLeaves = 0;
        float antisymmetry = 0.0f;
        ArrayList<Node> nodes = getRoot().depthFirstTraversal();
        nodeCount = nodes.size();
        for (Node node: nodes) {
            Node n = node;
            int depth = n.getDepth();
            nChildren += n.children.size();
            if (depth > maxDepth) {
                maxDepth = depth;
            }
            if (node.children.size() == 0) {
                // ie if this is a leaf
                nLeaves++;
                depthSum += depth;
                if (depth < minDepth) {
                    minDepth = depth;
                }
            } else {
                float epsilon = 0.0001f;
                float middleChild = node.children.size() / 2.0f - epsilon;
                int leftSize = 0, rightSize = 0;
                for (int i = 0; i < middleChild; i++) {
                    leftSize += node.children.get(i).getSubtreeSize();
                }
                for (int i = node.children.size() - 1; i > middleChild; i--) {
                    rightSize += node.children.get(i).getSubtreeSize();
                }
                antisymmetry += Math.abs(leftSize - rightSize);
            }
        }
        meanDepth = depthSum / (float) nLeaves;
        meanFanout = nChildren / (float) getNodeCount();
        symmetry = 1.0f - (antisymmetry / nodeCount);
    }

    public int getNodeCount() {
        if (nodeCount < 0) {
            calculateFeatures();
        }
        return nodeCount;
    }

    public int getMaxDepth() {
        if (maxDepth < 0) {
            calculateFeatures();
        }
        return maxDepth;
    }

    public int getMinDepth() {
        if (minDepth > 999) {
            calculateFeatures();
        }
        return minDepth;
    }

    public float getMeanDepth() {
        if (meanDepth < 0) {
            calculateFeatures();
        }
        return meanDepth;
    }

    public float getMeanFanout() {
        if (meanFanout < 0) {
            calculateFeatures();
        }
        return meanFanout;
    }

    public float getSymmetry() {
        // careful: symmetry could go below 0, so check something else
        // as the flag for whether features have already been calculated
        if (nodeCount < 0) {
            calculateFeatures();
        }
        return symmetry;
    }

    private void parseString(Node parent, StringTokenizer st) {

        while (st.hasMoreTokens()) {
            String currTok = st.nextToken().trim();

            if (currTok.equals("")) {
                // Tokenizer gave us an empty token, do nothing.

            } else if (currTok.equals("(") || currTok.equals("[")) {
                // The next token is the parent of a new subtree
                currTok = st.nextToken().trim();
                Node newNode = new Node(parent, currID++, currTok);
                parent.addChild(newNode);
                parseString(newNode, st);

            } else if (currTok.equals(")") || currTok.equals("]")) {
                // Finish this subtree
                return;

            } else {
                // An ordinary child node: add it to parent and continue.
                Node newNode = new Node(parent, currID++, currTok);
                parent.addChild(newNode);
            }
        }
        return;
    }

    public String toString() {
        return getRoot().toStringAsTree();
    }


    public double eval() {
        double res = 0.0;
        try {
            // do a string replace: / to p (meaning pdiv)
            if (program == null) {
                program = toString();
                program = program.replace('/', 'p');
            }
            res = (Double) js.eval(program);
        } catch (Exception ex) {
            ex.printStackTrace();
            System.err.println("woops! ");
        }
        return res;
    }

    public double fitness() {
        // use p here instead of /
        String target = "(define target (lambda (x y) (* x (+ x y))))";
        String candidate = toString();
        double fitness = 0.0;
        try {
            js.eval(target);
            for (int i = 0; i < 10; i++) {
                String setup = "(define x " + i / 10.0 + ")";
                js.eval(setup);
                for (int j = 0; j < 10; j++) {
                    String setup2 = "(define y " + j / 10.0 + ")";
                    js.eval(setup2);
                    fitness += abs(((Double) js.eval("(target x y)")) - eval());
                }
            }
            fitness /= 100.0;
        } catch (Exception ex) {
            ex.printStackTrace();
            System.err.println("woops! ");
        }
        return fitness;
    }


    public static void main(String args[]) {
        Tree trees[] = {
            new Tree("(* x y)"),
            new Tree("(+ (* x x) (* y y))"),
            new Tree("(* x (/ x y))")
        };
        for (Tree t: trees) {
            System.out.println(t);
            FeatureVector fv = new FeatureVector(t);
            System.out.println(fv);
            // t.eval();
            System.out.println("fitness: " + t.fitness());
        }
        return;
    }
}
