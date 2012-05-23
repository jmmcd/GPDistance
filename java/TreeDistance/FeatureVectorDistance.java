/*
 * FeatureVectorDistance.java
 *
 * Created November 3 2010
 *
 */
package TreeDistance;
import static java.lang.Math.*;

/**
 * Calculate many features for each tree, and compute the distance
 * between the feature vectors. Take the mean difference between
 * values, normalised as necessary, except for hashCode: here use the
 * discrete metric.
 *
 * @author James McDermott
 */
public class FeatureVectorDistance {

    int maxNodes;
    int maxDepth;
    int maxArity;

    /**
     * Construct a Feature Vector Distance object.
     *
     * @param _maxDepth maximum depth of the language.
     * @param _maxArity maximum arity of any node in the language.
     */
    public FeatureVectorDistance(int _maxDepth, int _maxArity) {
        // +1 because we count the root as depth 0, but it's "really" 1.
        maxDepth = _maxDepth + 1;
        maxArity = _maxArity;

        maxNodes = 0;
        for (int i = 0; i < maxDepth; i++) {
            maxNodes += pow(maxArity, i);
        }
    }

    /**
     * The main entry point to this class. Return a double representing
     * the distance between input trees.
     *
     * @param t
     * @param s a pair of Trees, the distance between them is to be
     * calculated.
     */
    public double fvd(Tree t, Tree s) {
        FeatureVector tfv = new FeatureVector(t), sfv = new FeatureVector(s);

        double distance =
            getNodeCount(tfv, sfv) +
            getMinDepth(tfv, sfv) +
            getMeanDepth(tfv, sfv) +
            getMaxDepth(tfv, sfv) +
            getSymmetry(tfv, sfv) +
            getMeanFanout(tfv, sfv) +
            getDiscreteMetricOnString(tfv, sfv);

		// System.out.println("FVD values:");
		// System.out.println(getNodeCount(tfv, sfv) + " " +
		// 				   getMinDepth(tfv, sfv) + " " +
		// 				   getMeanDepth(tfv, sfv) + " " +
		// 				   getMaxDepth(tfv, sfv) + " " +
		// 				   getSymmetry(tfv, sfv) + " " +
		// 				   getMeanFanout(tfv, sfv) + " " +
		// 				   getDiscreteMetricOnString(tfv, sfv)
		// 				   );
        distance /= 7.0;
        return distance;
    }

    public double getNodeCount(FeatureVector tfv, FeatureVector sfv) {
        return abs(tfv.nodeCount - sfv.nodeCount) / (float) maxNodes;
    }

    public double getMaxDepth(FeatureVector tfv, FeatureVector sfv) {
        return abs(tfv.depth - sfv.depth) / (float) maxDepth;
    }

    public double getMinDepth(FeatureVector tfv, FeatureVector sfv) {
        return abs(tfv.minDepth - sfv.minDepth) / (float) maxDepth;
    }

    public double getMeanDepth(FeatureVector tfv, FeatureVector sfv) {
        return abs(tfv.meanDepth - sfv.meanDepth) / (float) maxDepth;
    }

    public double getSymmetry(FeatureVector tfv, FeatureVector sfv) {
        return abs(tfv.symmetry - sfv.symmetry);
    }

    public double getMeanFanout(FeatureVector tfv, FeatureVector sfv) {
        return abs(tfv.meanFanout - sfv.meanFanout) / maxArity;
    }

    public double getDiscreteMetricOnString(FeatureVector tfv,
											FeatureVector sfv) {
        if (!tfv.str.equals(sfv.str)) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
}
