package TreeDistance.TreeEditDistance;

import java.util.List;
import java.util.ArrayList;

/* This is a basic operations cost set for the Zhang and Shasha
 * algorithm. All operations are 1.  If node labels match, the cost is
 * 0.
 *
 * INSERT-LICENCE-INFO
 */
public class OpsZhangShasha {

    public static final int INSERT = 0;
    public static final int DELETE = 1;
    public static final int RENAME = 2;
    public static final int IDENTITY = 3;

    private ArrayList<TreeEditOperation> operations
        = new ArrayList<TreeEditOperation>();

    public OpsZhangShasha() {
        operations.add(INSERT, new BasicInsert());
        operations.add(DELETE, new BasicDelete());
        operations.add(RENAME, new BasicRename());
    }

    /** Returns the nth operation.  For convenience, the operation IDs are
     *  provided as named int constants:
     */
    public TreeEditOperation getOp(int opID) {
        return operations.get(opID);
    }

    public List<TreeEditOperation> getOperations() {
        return operations;
    }
}
