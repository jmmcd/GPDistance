package TreeDistance.TreeEditDistance;

import java.util.List;
import java.util.ArrayList;

/** This class stores the sequence of operations that transform one
 * tree into another.  The two trees are stored in this data
 * structure.
 *
 * INSERT-LICENCE-INFO
 */
public class Transformation {

    private TreeDefinition aTree;
    private TreeDefinition bTree;
    private double totalCost;
    private ArrayList<TreeEditOperation> operationsList
        = new ArrayList<TreeEditOperation>();

    /** Returns the total cost of the transformation */
    public double getCost() {
        return totalCost;
    }

    /** Returns the total cost of the transformation */
    public void setCost(double _cost) {
        totalCost = _cost;
    }

    /** Returns the (least costly) list of operation required to
     * transform TreeA into TreeB.
     */
    public List<TreeEditOperation> getOperations() {
        return operationsList;
    }
}
