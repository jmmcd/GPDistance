package TreeDistance.TreeEditDistance;

/* Implements a basic insert operation with cost 1.
 *
 * INSERT-LICENCE-INFO
 */
public class BasicInsert extends TreeEditOperation {

    public BasicInsert() {
	super.opName = "INSERT";
    }

    public double getCost(int aNodeID, int bNodeID,
			  TreeDefinition aTree,
			  TreeDefinition bTree) {
	return 1;
    }
}
