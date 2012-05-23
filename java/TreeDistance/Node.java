/*
 * Node.java
 *
 * Created November 3 2010
 *
 */
package TreeDistance;

import java.util.ArrayList;

public class Node {
    public Node parent;
    public int id;
    public String label;
    public ArrayList<Node> children;
    public int subtreeSize;
    private int depth;

    public Node(Node _parent,
                int _id,
                String _label) {
        parent = _parent;
        id = _id;
        label = _label;
        children = new ArrayList<Node>();
        subtreeSize = -1;
        depth = -1;
    }

    public void addChild(Node child) {
        children.add(child);
    }

    public String toString() {
        return label;
    }

    public Tree cloneAsTree() {
        return new Tree(toStringAsTree());
    }

    public int getDepth() {
        if (depth == -1) {
            // Travel back to root to calculate depth
            int _depth = 0;
            Node n = this;
            while (!n.parent.label.equals("holder")) {
                n = n.parent;
                _depth++;
            }
            depth = _depth;
        }

        return depth;
    }


    public String toStringAsTree() {
        if (children.size() > 0) {
            String retval = "(" + this;
            for (Node child: children) {
                retval += " " + child.toStringAsTree();
            }
            retval += ")";
            return retval;
        } else {
            return this.toString();
        }
    }

    /**
     * Do a depth-first traversal of the tree starting at a given node.
     * @return a list of Nodes in depth-first order.
     */
    public ArrayList<Node> depthFirstTraversal() {
        ArrayList<Node> retval = new ArrayList<Node>();
        retval.add(this);
        for (Node child: children) {
            retval.addAll(child.depthFirstTraversal());
        }
        return retval;
    }

    public int getSubtreeSize() {
        if (subtreeSize == -1) {
            subtreeSize = 1;
            for (Node child: children) {
                subtreeSize += child.getSubtreeSize();
            }
        }
        return subtreeSize;
    }

    public ArrayList<Node> getAncestors() {
        Node n = this;
        ArrayList<Node> retval = new ArrayList<Node>();
        retval.add(n);
        while (!n.parent.label.equals("holder")) {
            n = n.parent;
            retval.add(n);
        }
        return retval;
    }


}
