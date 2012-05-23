package TreeDistance.TreeEditDistance;

/* This is an abstract definition of what is needed to define a Tree
 * Edit Operation.  You can create subclasses of this to meet your
 * application needs if more complicated inferences on edit costs are
 * required.
 *
 * INSERT-LICENCE-INFO
 */
public abstract class TreeEditOperation {

    String opName = "";

    public String getName() {
	return opName;
    }

    public abstract double getCost(int aNodeID, int bNodeID,
				   TreeDefinition aTree,
				   TreeDefinition bTree);

}
