package TreeDistance.TreeEditDistance;

import java.util.Hashtable;
import java.util.ArrayList;

/*
 * INSERT-LICENCE-INFO
 */
public class TestZhangShasha {

    public static void main
		(String argv []) throws java.io.IOException  {

		if (argv.length != 2) {
			System.out.println("Usage TestZhangShasha <tree1> <tree2>");
			return;
		}

		TreeDefinition aTree = CreateTreeHelper.makeTree(argv[0]);
		System.out.println("The tree is: \n"+aTree);
		TreeDefinition bTree = CreateTreeHelper.makeTree(argv[1]);
		System.out.println("The tree is: \n"+bTree);

		ComparisonZhangShasha treeCorrector = new ComparisonZhangShasha();
		OpsZhangShasha costs = new OpsZhangShasha();
		Transformation transform = treeCorrector.findDistance(aTree, bTree, costs);
		System.out.println("Distance: "+transform.getCost());
		System.out.println("Trace: " + transform.getOperations());
    }

		// Hashtable<Character, Integer> arities = new Hashtable<Character, Integer>();
		// arities.put('+', 2);
		// arities.put('*', 2);
		// arities.put('x', 0);
		// arities.put('1', 0);
        // char [] input1 = {'+', '*', 'x', 'x', '1'};
        // char [] input2 = {'+', '*', 'x', '1', '1'};

		// CreateTreeHelper cth = new CreateTreeHelper();
		// Tree t1 = cth.makeTreeFromPreorder(input1, arities);
		// Tree t2 = cth.makeTreeFromPreorder(input2, arities);
		// System.out.println(t1);
		// System.out.println(t2);

}
