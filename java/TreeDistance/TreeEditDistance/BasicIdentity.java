package TreeDistance.TreeEditDistance;

/* Implements a basic identity operation with cost 0.
 *
 * INSERT-LICENCE-INFO
 */
public class BasicIdentity extends TreeEditOperation {

    public BasicIdentity() {
	super.opName = "IDENTITY";
    }

    public double getCost(int aNodeID, int bNodeID,
			  TreeDefinition aTree,
			  TreeDefinition bTree) {
	return 0;
    }
}
