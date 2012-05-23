package TreeDistance.TreeEditDistance;

/* Implements a basic delete operation with cost 1.
 *
 * INSERT-LICENCE-INFO
 */
public class BasicDelete extends TreeEditOperation {

    public BasicDelete() {
	super.opName = "DELETE";
    }

    public double getCost(int aNodeID, int bNodeID,
			  TreeDefinition aTree,
			  TreeDefinition bTree) {
	return 1;
    }
}
