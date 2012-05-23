package TreeDistance.TreeEditDistance;

import java.util.Hashtable;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.TreeSet;


/** This defines necessary operations on tree required for the tree
 * algorithms in this package.  Subclasses handle the particular
 * representation of the contents of the tree. This parent superclass
 * looks after generic tree functionality like ordering nodes.
 *
 * INSERT-LICENCE-INFO
 */
public abstract class TreeDefinition {

    //Just two constants for defining ordering.
    public static final int POSTORDER = 0;
    public static final int PREORDER = 1;
    private int chosenOrder = 0; //default POSTORDER

    //root node label
    private String root = "";

    private Hashtable<Integer, String> iDsToLabel
        = new Hashtable<Integer, String>() ;
    private Hashtable<String, Integer> labelToIDs
        = new Hashtable<String, Integer>() ;
    private Hashtable<Integer, ArrayList<Integer>> treeStructureIDs
        = new Hashtable<Integer, ArrayList<Integer>>();

    public void setRoot(String _root) {
        root = _root;
    }

    public String getRoot() {
        return root;
    }

    /** Returns the NodeID of the root */
    public int getRootID() {
        return labelToIDs.get(root);
    }

    /** This is the ordering used to number the nodes */
    public void setOrder(int order) {
        chosenOrder = order;
    }

    /** This is the ordering used to number the nodes */
    public int getOrder() {
        return chosenOrder;
    }

    public void orderNodes(int ordering) {
        //Set ordering
        if (ordering == POSTORDER) {
            setPostOrdering(0, root);
            setOrder(POSTORDER);
        }
        else { //PREORDER
            //      setPreOrdering(0, root);
            setOrder(PREORDER);
        }

        //  System.out.println("TreeDefinition.init: PostOrder IDs");
        //  for (String label : new TreeSet<String>(labelToIDs.keySet())) {
        //      System.out.println("TreeDefinition.init: Label: "+label+", id: "+
        //                 labelToIDs.get(label));
        //  }

        //
        //Create version of tree that just uses index numbers
        //
        for (String parent : getNodes()) {

            //      System.out.println("TreeDefinition: parent: "+parent);

            ArrayList<Integer> indexedChildren
                = new ArrayList<Integer>();
            for (String child : getChildren(parent)) {

                //      System.out.println("TreeDefinition: child: "+child);

                indexedChildren.add(getNodeID(child));
            }
            treeStructureIDs.put(getNodeID(parent),
                                 indexedChildren);
        }
    }

    public int setPostOrdering(int counter, String aNodeLabel) {
        int internalCounter = counter;

        //examine children
        for (String child: getChildren(aNodeLabel)) {
            internalCounter = setPostOrdering(internalCounter, child);
        }

        //set new nodeID for this node (set to counter+1)
        putLabel(internalCounter+1, aNodeLabel);
        putNodeID(aNodeLabel, internalCounter+1);

        return internalCounter+1;
    }


    /** This is provides the string label for the node we're matching.
     * In some cases, the value of the string may depend on properties
     * of the node in addition to the actual node label. */
    public String getLabelForMatching(int nodeID) {
        return getLabel(nodeID);
    }

    public String getLabel(int nodeID) {
        return iDsToLabel.get(nodeID);
    }

    public int getNodeID(String nodeLabel) {
        return labelToIDs.get(nodeLabel);
    }

    public void putLabel(int nodeID, String nodeLabel) {
        iDsToLabel.put(nodeID, nodeLabel);
    }

    public void putNodeID(String nodeLabel, int nodeID) {
        labelToIDs.put(nodeLabel, nodeID);
    }

    /* returns a collection of nodes referred to by their original
     * (unique) nodel label.
     */
    public abstract Collection<String> getNodes();

    /* returns a collection of children referred to by their original
     * (unique) nodel label.
     */
    public abstract List<String> getChildren(String aNodeLabel);

    /** Returns the children of the node given as a parameter. */
    public Collection<Integer> getChildrenIDs(int nodeID) {
        return treeStructureIDs.get(nodeID);
    }
    public boolean isLeaf(int nodeID) {
        return (getChildrenIDs(nodeID).size() == 0);
    }

    public int getNodeCount() {
        return treeStructureIDs.keySet().size();
    }

    public abstract String toString();
}
