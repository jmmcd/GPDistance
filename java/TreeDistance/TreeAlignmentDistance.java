/*
 * TreeAlignmentDistance.java
 *
 * Created on 29 September 2009
 *
 * Implements the distance proposed by Shin-Yee Lu, or Stanley M.
 * Selkow, or Kuo-Chung Tai, and cited by Ekart and Nemeth and then
 * used by Vanneschi et al.
 *
 */
package TreeDistance;

import java.util.HashMap;
import java.util.List;

/**
 *
 * @author James McDermott
 */
public class TreeAlignmentDistance {

    boolean useDiscrete;
    HashMap<String, Integer> codes;
    int z = 2;
    double weight;
	double k = 1.0;

    // This constructor will use the discrete metric on the
    // primitives: if they differ, then the distance is 1.
    public TreeAlignmentDistance() {
        useDiscrete = true;
        weight = 0.5;
    }

    // This constructor will use the discrete metric on the
    // primitives: if they differ, then the distance is 1. Specify the
    // depth-weighting parameter.
    public TreeAlignmentDistance(double _weight) {
        useDiscrete = true;
        weight = _weight;
    }

    // This constructor allows you to specify a coding function as a
    // hashmap over primitives. Should not map to zero since that's
    // reserved for null. Specify depth-weighting parameter also.
    public TreeAlignmentDistance(HashMap<String, Integer> _codes, double _weight) {
        useDiscrete = false;
        codes = _codes;
        weight = _weight;
    }

	// This constructor will use a coding function which starts at 1
	// for the first primitive and increments for each subsequent one.
	// This was assumed in the Vanneschi work, where for each
	// primitive, code = arity + 1. Code zero is reserved for null.
	// Specify depth-weighting parameter also.
	public TreeAlignmentDistance(List<String> primitives, double _weight) {
		useDiscrete = false;
		codes = new HashMap<String, Integer>();
		for (int i = 0; i < primitives.size(); i++) {
			codes.put(primitives.get(i), i + 1);
		}
        weight = _weight;
	}

	// // Public method for distance between individuals. We assume that 
	// // we only want the distance between individuals' first trees.
	// public double getDistance(GPIndividual i1, GPIndividual i2) {
	// 	return getDistance(i1.trees[0], i2.trees[0]);
	// }

    // The public method for calculating distance between trees.
    public double getDistance(Tree t1, Tree t2) {
        return distance(t1.getRoot(), t2.getRoot(), k);
    }

    // This recursively calculates the distance.
    private double distance(Node t1, Node t2, double k) {

        double retval = 0.0;
        int m;

        // Calculate the first part: distance between primitives
        if (t1 == null && t2 == null) {
            return retval;
        } else if (t1 == null) {
            retval += d(null, t2.toString());
            m = t2.children.size();
        } else if (t2 == null) {
            retval += d(t1.toString(), null);
            m = t1.children.size();
        } else {
            retval += d(t1.toString(), t2.toString());
            m = Math.max(t1.children.size(), t2.children.size());
        }

        // Sum the second part: distance between each pair of children
        for (int i = 0; i < m; i++) {
            Node ptr1 = null;
            Node ptr2 = null;
            if (t1 != null && i < t1.children.size()) {
                ptr1 = t1.children.get(i);
            }
            if (t2 != null && i < t2.children.size()) {
                ptr2 = t2.children.get(i);
            }
            // Recursive call.
            retval += k * distance(ptr1, ptr2, k * weight);
        }
        return retval;
    }

    // Distance between a pair of primitives.
    private double d(String e1, String e2) {
        if (useDiscrete == false) {
            // Coding function
            return Math.pow(Math.abs(code(e1) - code(e2)), z);
        } else {
            // Discrete metric on primitives
            if (e1 == null && e2 == null) {
                return 0;
            } else if (e1 == null || e2 == null) {
                return 1;
            } else if (e1.equals(e2)) {
                return 0;
            } else {
                return 1;
            }
        }
    }

    private int code(String e) {
		if (e == null) {
			return 0;
		}
        return codes.get(e);
    }
}
