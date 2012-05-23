package TreeDistance.TreeEditDistance;

import java.util.Hashtable;
import java.util.ArrayList;

/*
 * INSERT-LICENCE-INFO
 */
public class TestTreeDefinition {

    public static BasicTree aBasicTree;


    public static void main
	(String argv []) throws java.io.IOException  {

	System.out.println("Testing Tree Definition.");

	//A Tree
	Hashtable<String, ArrayList<String>> testTree
	    = new Hashtable<String, ArrayList<String>>();

	//a-[b,c]
	ArrayList<String> aChildren =
	    new ArrayList<String>();
	aChildren.add("b");
	aChildren.add("c");
	testTree.put("a", aChildren);

	testTree.put("b", new ArrayList<String>());

	//c-[d, e, f]
	ArrayList<String> cChildren =
	    new ArrayList<String>();
	cChildren.add("d");
	cChildren.add("e");
	cChildren.add("f");
	testTree.put("c", cChildren);

	testTree.put("d", new ArrayList<String>());
	testTree.put("e", new ArrayList<String>());
	testTree.put("f", new ArrayList<String>());

	aBasicTree =
	    new BasicTree(testTree, "a", BasicTree.POSTORDER);

	//Test Output
	System.out.println("Static Test Tree: \n");
	System.out.println("The number of nodes are: "+
			   aBasicTree.getNodeCount());
	System.out.println("The tree is: \n"+aBasicTree);


	//Now test CreateTreeHelper
	TreeDefinition bBasicTree = null;

	if (argv.length == 1) {
	    bBasicTree = CreateTreeHelper.makeTree(argv[0]);
	    System.out.println("Input Tree: \n");
	    System.out.println("The number of nodes are: "+
			       bBasicTree.getNodeCount());
	    System.out.println("The tree is: \n"+bBasicTree);
	}
    }

}
