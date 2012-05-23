package TreeDistance.TreeEditDistance;

import TreeDistance.*;

// Copyright (C) 2004 Stephen Wan:
// http://web.science.mq.edu.au/~swan/howtos/treedistance/package.html

import java.util.Hashtable;
import java.util.ArrayList;
import java.util.Stack;

/*
 * INSERT-LICENCE-INFO
 */
public class CreateTreeHelper {

    private int tmpNodeID;

    public CreateTreeHelper() {
    }
    /* This takes a String describing a tree and converts it into a
     * TreeDefinition.  The format of the string is a series of edges
     * represented as pairs of string separated by semi-colons.  Each
     * pair is comma separated.  The first substring in the pair is
     * the parent, the second is the child.  The first edge parent
     * must be the root of the tree.
     *
     * For example: "a-b;a-c;c-d;c-e;c-f;"
     */
    public static TreeDefinition makeTree(String treeSpec)  {
		return makeTree(treeSpec, null);
    }


    public class ZSNode {
		public int id, arityRemaining;
		public String uniqueLabel;
        public String label;
		public ZSNode(int _id, String _label, int _arity) {
			id = _id;
			label = _label;
			arityRemaining = _arity;
			uniqueLabel = label + ":" + id;
		}
    }

    public BasicTree makeTreeFromPreorder(String[] input, Hashtable<String, Integer> arities) {

		Hashtable<String, ArrayList<String>> aTree
			= new Hashtable<String, ArrayList<String>>();

        Stack<ZSNode> parents = new Stack<ZSNode>();

		ZSNode currParent = null;
		String rootLabel = null;

		int id = 0;
        for (String currTok: input) {
			int arity = arities.get(currTok);
			ZSNode currNode = new ZSNode(id++, currTok, arity);

			// first token: it's the root
			if (parents.empty()) {
				// FIXME hack hack! We add this node in case it has no children
				// so addEdge() will never be triggered.
				addNode(currNode.uniqueLabel, aTree);
				rootLabel = currNode.uniqueLabel;
			} else {
				// not first token: add incoming edge
				addEdge(currParent.uniqueLabel, currNode.uniqueLabel, aTree);
				currParent.arityRemaining -= 1;
			}

			// this node is a parent: save previous parent and use this
			if (arity > 0) {
				parents.push(currParent);
				currParent = currNode;
			}

			// current parent has sufficient children: restore previous.
			// if currParent == null, this is the first node in the tree,
			// and it's a terminal
			if (currParent != null && currParent.arityRemaining == 0) {
				currParent = parents.pop();
			}
		}
		BasicTree aBasicTree =
			new BasicTree(aTree, rootLabel, BasicTree.POSTORDER);
		return aBasicTree;
    }

    public BasicTree makeBasicTreeFromTree(Tree t) {
		Hashtable<String, ArrayList<String>> aTree
			= new Hashtable<String, ArrayList<String>>();

        for (Node node: t.getRoot().depthFirstTraversal()) {
            addNode(node.label + ":" + node.id, aTree);
            if (!node.parent.label.equals("holder")) {
                addEdge(node.parent.label + ":" + node.parent.id,
                        node.label + ":" + node.id,
                        aTree);
            }
        }

        Node parent = t.getRoot();
        BasicTree aBasicTree =
			new BasicTree(aTree, parent.label + ":" + parent.id, BasicTree.POSTORDER);

        return aBasicTree;
    }

    // public BasicTree makeTreeFromECJGPTree(GPTree gpt) {
    //     tmpNodeID = 0;

    //     Hashtable<String, ArrayList<String>> aTree
    //         = new Hashtable<String, ArrayList<String>>();

    //     // use dummy arity 0, because it's not needed here.
    //     ZSNode currNode = new ZSNode(tmpNodeID++, gpt.child.toString(), 0);
    //     addNode(currNode.uniqueLabel, aTree);
    //     addEdges(aTree, gpt.child, currNode);

    //     String rootLabel = gpt.child.toString() + Integer.toString(gpt.getRoot().getID());
    //     BasicTree aBasicTree =
	// 		new BasicTree(aTree, currNode.uniqueLabel, BasicTree.POSTORDER);

    //     return aBasicTree;
    // }


    // // Recursive function for adding edges.
    // private void addEdges(Hashtable<String, ArrayList<String>> aTree, ZSNode n,
    //                       ZSNode parent) {
    //     // Iterate through children
    //     for (String child: aTree.get(n)) {
	// 		ZSNode currNode = new ZSNode(tmpNodeID++, child, 0);
    //         addNode(currNode.uniqueLabel, aTree);
    //         addEdge(parent.uniqueLabel, currNode.uniqueLabel, aTree);
    //         // Recursive call
    //         addEdges(aTree, child, currNode);
    //     }
    // }


    public void printEdge(String s, String t) {
		System.out.println(s + " - " + t);
    }

    /* This takes a String describing a tree and converts it into a
     * TreeDefinition.  The format of the string is a series of edges
     * represented as pairs of string separated by semi-colons.  Each
     * pair is comma separated.  The first substring in the pair is
     * the parent, the second is the child.
     *
     * For example: "a-b;a-c;c-d;c-e;c-f;"
     */
    public static TreeDefinition makeTree(String treeSpec, String rootID)  {

		//A Tree
		Hashtable<String, ArrayList<String>> aTree
			= new Hashtable<String, ArrayList<String>>();

		String root = rootID;

		String[] edges = treeSpec.split(";");
		for (String edge: java.util.Arrays.asList(edges)) {

			System.out.println("CreateTreeHelper: Examining edge: "+edge);

			String[] nodes = edge.split("-");
			if (nodes.length > 1) {
				addEdge(nodes[0], nodes[1], aTree);
			} else {
				addNode(nodes[0], aTree);
			}
			if (root == null) {
				root = nodes[0];
			}
		}

		BasicTree aBasicTree =
			new BasicTree(aTree, root, BasicTree.POSTORDER);

		return aBasicTree;
    }

    /** This adds a node if necessary to the tree definition.
     */
    protected static void addNode(String nodeLabel,
								  Hashtable<String, ArrayList<String>> treeStructure) {
		if (!treeStructure.containsKey(nodeLabel)) {
			treeStructure.put(nodeLabel, new ArrayList<String>());
		}
    }

    /** This adds the edge (and nodes if necessary) to the tree
     * definition .
     */
    protected static void addEdge(String parentLabel, String childLabel,
								  Hashtable<String, ArrayList<String>> treeStructure) {
		addNode(parentLabel, treeStructure);

		treeStructure.get(parentLabel).add(childLabel);

		addNode(childLabel, treeStructure);
    }

    public static void testSExpressions() {
		Hashtable<String, Integer> arities = new Hashtable<String, Integer>();
		arities.put("+", 2);
		arities.put("*", 2);
		arities.put("x", 0);
		arities.put("1", 0);
        String [] input = {"+", "*", "1", "1", "*", "x", "1"};

        System.out.println("Testing S-Expression and Tree Code.");
        CreateTreeHelper cth = new CreateTreeHelper();
        Tree nt1 = new Tree("(+ (* 1 1) (* x x))");
        Tree nt2 = new Tree("(+ x (* x x))");
        TreeDefinition t1 = cth.makeBasicTreeFromTree(nt1);
        TreeDefinition t2 = cth.makeBasicTreeFromTree(nt2);
		TreeDefinition t3 = cth.makeTreeFromPreorder(input, arities);

        System.out.println(t1);
        System.out.println(t2);
        System.out.println(t3);

		ComparisonZhangShasha treeCorrector = new ComparisonZhangShasha();
		OpsZhangShasha costs = new OpsZhangShasha();
		Transformation transform = treeCorrector.findDistance(t1, t2, costs);
		System.out.println("Distance: "+transform.getCost());

        transform = treeCorrector.findDistance(t1, t3, costs);
		System.out.println("Distance: "+transform.getCost());
    }

    public static void main(String[] args) {
		Hashtable<String, Integer> arities = new Hashtable<String, Integer>();
		arities.put("+", 2);
		arities.put("*", 2);
		arities.put("x", 0);
		arities.put("1", 0);
        String [] input1 = {"+", "*", "x", "x", "1"};
        String [] input2 = {"+", "*", "x", "1", "1"};

		CreateTreeHelper cth = new CreateTreeHelper();
		BasicTree t1 = cth.makeTreeFromPreorder(input1, arities);
		BasicTree t2 = cth.makeTreeFromPreorder(input2, arities);
		System.out.println(t1);
		System.out.println(t2);

        testSExpressions();
    }
}
