package TreeDistance.TreeEditDistance;

import java.util.Hashtable;
import java.util.List;
import java.util.ArrayList;
import java.util.TreeSet;

/** This is an implementation of the Zhang and Shasha algorithm as
 * described in [FIXME]
 *
 * SWAN 2007-11-01: I'm pretty sure this code comes from:
 * http://www.cs.queensu.ca/TechReports/Reports/1995-372.pdf and
 *"http://www.inf.unibz.it/dis/teaching/ATA/ata7-handout-1x1.pdf"
 * INSERT-LICENCE-INFO
 */
public class ComparisonZhangShasha {

    //"Dynamic Programming" Table.
    //use function setFD to access it.
    //Each call to findDistance will change these tables.  But each
    //call is independent (and reinitialises this) so the side effect
    //has no real consequence.  ie.  There are NO public side effects.
    Hashtable<Integer, Hashtable<Integer, Double>> fD = null;
    private double[][] distance = null;

    public Transformation findDistance
        (TreeDefinition aTree,
         TreeDefinition bTree,
         OpsZhangShasha ops) {

        //This is initialised to be n+1 * m+1.  It should really be n*m
        //but because of java's zero indexing, the for loops would
        //look much more readable if the matrix is extended by one
        //column and row.  So, distance[0,*] and distance[*,0] should
        //be permanently zero.
        distance =
            new double[aTree.getNodeCount()+1][bTree.getNodeCount()+1];

        //Preliminaries
        //1. Find left-most leaf and key roots
        Hashtable<Integer, Integer> aLeftLeaf = new Hashtable<Integer, Integer>();
        Hashtable<Integer, Integer> bLeftLeaf = new Hashtable<Integer, Integer>();
        ArrayList<Integer> aTreeKeyRoots = new ArrayList<Integer>();
        ArrayList<Integer> bTreeKeyRoots = new ArrayList<Integer>();

        findHelperTables(aTree, aLeftLeaf, aTreeKeyRoots, aTree.getRootID());

        findHelperTables(bTree, bLeftLeaf, bTreeKeyRoots, bTree.getRootID());

        // System.out.println("aTree Keyroots");
        for (Integer aKeyroot : aTreeKeyRoots) {
            // System.out.println(aKeyroot);
        }

        // System.out.println("l(aTree)");
        for (Integer key : new TreeSet<Integer>(aLeftLeaf.keySet())) {
            // System.out.println(key+": "+aLeftLeaf.get(key));
        }

        // System.out.println("bTree Keyroots");
        for (Integer bKeyroot : bTreeKeyRoots) {
            // System.out.println(bKeyroot);
        }

        // System.out.println("l(bTree)");
        for (Integer key : new TreeSet<Integer>(bLeftLeaf.keySet())) {
            // System.out.println(key+": "+bLeftLeaf.get(key));
        }

        //Comparison
        for (Integer aKeyroot : aTreeKeyRoots) { //aKeyroot loop
            // System.out.println("aKeyroot: "+aKeyroot);

            for (Integer bKeyroot : bTreeKeyRoots) { //bKeyroot loop
                // System.out.println("bKeyroot: "+bKeyroot);


                // Reinitialise forest-distance dynamic programming table.
                fD = new Hashtable<Integer, Hashtable<Integer, Double>>();

                setFD(aLeftLeaf.get(aKeyroot),bLeftLeaf.get(bKeyroot),0.0d);

                //for all descendents of aKeyroot: i
                for (int i=aLeftLeaf.get(aKeyroot);i <= aKeyroot;i++) {
                    setFD(i,
                          bLeftLeaf.get(bKeyroot)-1,
                          getFD(i-1,bLeftLeaf.get(bKeyroot)-1)+
                          ops.getOp(OpsZhangShasha.DELETE).getCost(i,0,aTree,bTree));

                    //          //trace
                    // seeFD(i,bLeftLeaf.get(bKeyroot)-1);
                }

                //for all descendents of bKeyroot: j
                for (int j=bLeftLeaf.get(bKeyroot);j <= bKeyroot;j++) {

                    setFD(aLeftLeaf.get(aKeyroot)-1,j,
                          getFD(aLeftLeaf.get(aKeyroot)-1,j-1)+
                          ops.getOp(OpsZhangShasha.INSERT).getCost(0,j,aTree,bTree));

                    //          //trace
                    // seeFD(aLeftLeaf.get(aKeyroot)-1,j);
                }

                // System.out.println("BEFORE cross checks: ");
                // seeFD();

                //for all descendents of aKeyroot: i
                for (int i=aLeftLeaf.get(aKeyroot);i<=aKeyroot;i++) {
                    //              System.out.println("i: "+i);

                    //for all descendents of bKeyroot: j
                    for (int j=bLeftLeaf.get(bKeyroot);j<=bKeyroot;j++) {
                        //               System.out.println("j: "+j);

                        //Start Trace
                                  // System.out.println
                                      // ("DEL: "+
                                      //  (getFD(i-1,j)+
                                      //  ops.getOp(OpsZhangShasha.DELETE)
                                      //   .getCost(i,0,aTree,bTree)));

                                  // System.out.println
                                      // ("INS: "+
                                      //  (getFD(i,j-1)+
                                      //  ops.getOp(OpsZhangShasha.INSERT)
                                      //   .getCost(0,j,aTree,bTree)));

                        //End Trace
                        double min =  //This min compares del vs ins
                            java.lang.Math.min
                            (//Option 1: Delete node from aTree
                             getFD(i-1,j)+
                             ops.getOp(OpsZhangShasha.DELETE)
                             .getCost(i,0,aTree,bTree),

                             //Option 2: Insert node into bTree
                             getFD(i,j-1)+
                             ops.getOp(OpsZhangShasha.INSERT)
                             .getCost(0,j,aTree,bTree)
                             );

                        //           System.out.println("Min: "+min);

                        if ((aLeftLeaf.get(i) == aLeftLeaf.get(aKeyroot))
                            &&
                            (bLeftLeaf.get(j) == bLeftLeaf.get(bKeyroot)))
                            {

                                                 // System.out.println("This is a Left-branch node.");

                                             // System.out.println
                                                 // ("REN: "+
                                                 //  (getFD(i-1,j-1) +
                                                 //  ops.getOp(OpsZhangShasha.RENAME)
                                                 //   .getCost(i,j,aTree,bTree)));

                                distance[i][j] =
                                    java.lang.Math.min
                                    (min,
                                     getFD(i-1,j-1) +
                                     ops.getOp(OpsZhangShasha.RENAME)
                                     .getCost(i,j,aTree,bTree)
                                     );

                                // System.out.println("D["+i+"]["+j+"]:"+
                                // distance[i][j]);

                                setFD(i,j,distance[i][j]);

                                //trace
                                // seeFD(i,j);
                            }
                        else {

                            // System.out.println("Forest Situation");

                            // System.out.println
                            //     ("REN: "+
                            //      (getFD(aLeftLeaf.get(i)-1,
                            //             bLeftLeaf.get(j)-1)+
                            //       distance[i][j] ));

                            setFD(i,j,
                                  java.lang.Math.min
                                  (min,
                                   getFD(aLeftLeaf.get(i)-1,
                                         bLeftLeaf.get(j)-1)+
                                   distance[i][j]
                                   ));

                            // seeFD(i, j);
                        }
                    }
                }
                //trace
                // seeFD();
            }
        }

        //  //Return result
        //  for (int i=0;i<distance.length;i++) {
        //      for (int j=0;j<distance[i].length;j++) {
        //      System.out.print(distance[i][j]+" ");
        //      }
        //      System.out.println("");
        //  }

        //      System.out.println("Distance: "
        //             +distance[aTree.getNodeCount()][bTree.getNodeCount()]);

        Transformation transform = new Transformation();
        transform.setCost(distance[aTree.getNodeCount()][bTree.getNodeCount()]);
        return transform;

    }

    /** The initiating call should be to the root node of the tree.
     * It fills in an nxn (hash) table of the leftmost leaf for a
     * given node.  It also compiles an array of key roots. The
     * integer values IDs must come from the post-ordering of the
     * nodes in the tree.
     */
    private void findHelperTables
        (TreeDefinition someTree,
         Hashtable<Integer, Integer> leftmostLeaves,
         ArrayList<Integer> keyroots,
         int aNodeID) {

        findHelperTablesRecurse(someTree,leftmostLeaves,keyroots,aNodeID);

        //add root to keyroots
        keyroots.add(aNodeID);

        //add boundary nodes
        leftmostLeaves.put(0,0);

    }
    private void findHelperTablesRecurse
        (TreeDefinition someTree,
         Hashtable<Integer, Integer> leftmostLeaves,
         ArrayList<Integer> keyroots,
         int aNodeID) {

        //If this is a leaf, then it is the leftmost leaf
        if (someTree.isLeaf(aNodeID)) {
            leftmostLeaves.put(aNodeID, aNodeID);}
        else {
            boolean seenLeftmost = false;
            for (Integer child : someTree.getChildrenIDs(aNodeID)) {
                findHelperTablesRecurse(someTree, leftmostLeaves,
                                        keyroots, child);
                if (!seenLeftmost) {
                    leftmostLeaves.put(aNodeID,
                                       leftmostLeaves.get(child));
                    seenLeftmost = true;}
                else {
                    keyroots.add(child);
                }
            }
        }
    }

    /** Returns a String for trace writes */
    private void seeFD(int a, int b) {

        System.out.println("["+a+"],["+b+"]: "+getFD(a,b));
    }

    /** Returns a String for trace writes */
    private void seeFD() {
        System.out.println("Forest Distance");
        //Return result
        for (Integer i : new TreeSet<Integer>(fD.keySet())) {
            System.out.print(i+": ");
            for (Integer j : new TreeSet<Integer>(fD.get(i).keySet())) {
                System.out.print(fD.get(i).get(j)+"("+j+")  ");
            }
            System.out.println("");
        }
    }

    /** This returns the value in the cell of the ForestDistance table
     */
    private double getFD(int a, int b) {

        if (!fD.containsKey(a)) {
            //          System.out.println("getFD: creating new aStr entry.");
            fD.put(a, new Hashtable<Integer, Double>());
        }

        Hashtable<Integer, Double> row = fD.get(a);
        if (!row.containsKey(b)) {
            //      System.out.println("creating new bStr entry.");
            row.put(b, 0.0);
        }
        return row.get(b);
    }


    /** This sets the value in the cell of the ForestDistance table
     */
    private void setFD(int a, int b,
                       double aValue) {

        if (!fD.containsKey(a)) {
            //          System.out.println("setFD: creating new aStr entry.");
            fD.put(a, new Hashtable<Integer, Double>());
        }

        Hashtable<Integer, Double> row = fD.get(a);
        row.put(b, aValue);

        //      for (String key: fD.keySet()) {
        //          System.out.println("FD key: "+key);
        //      for (String key2: fD.get(key).keySet()) {
        //      System.out.println("FD key2: "+key2);
        //      }
        //      }
    }

}
