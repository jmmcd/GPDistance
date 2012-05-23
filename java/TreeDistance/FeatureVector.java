package TreeDistance;

public class FeatureVector {
    public int nodeCount;
    public int depth;
    public int minDepth;
    public float meanDepth;
    public float symmetry;
    public float meanFanout;

    public String str;

    public FeatureVector(Tree t) {
        nodeCount = t.getNodeCount();
        depth = t.getMaxDepth();
        minDepth = t.getMinDepth();
        meanDepth = t.getMeanDepth();
        symmetry = t.getSymmetry();
        meanFanout = t.getMeanFanout();
		str = t.toString();
    }

    public String toString() {
        return ("Nodecount = " + nodeCount + "; depth = " + depth + "; minDepth = " +
                minDepth + "; meanDepth = " + meanDepth + "; symmetry = " + symmetry +
                "; meanFanout = " + meanFanout);
    }
}

