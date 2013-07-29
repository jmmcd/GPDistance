/*
 * Language represents the terminal and function sets of the GP
 * system.
 *
 */

package GP;

import TreeDistance.*;

import java.util.*;
import static java.lang.Math.*;

public class Language {

    // Map from Strings (representing functions and terminals) to arities.
    public HashMap<String, Integer> arities;
    ArrayList<String> terminals;
    ArrayList<String> functions;

    // Number of functions and terminals.
    public int F, T;
    public int maxArity, maxDepth;

    // PT is the probability of choosing a terminal
    float PT;

    public static HashMap<String, Integer> defaultArities;
    static {
        defaultArities = new HashMap<String, Integer>();
        defaultArities.put("*", 2);
        defaultArities.put("+", 2);
        defaultArities.put("/", 2);
        defaultArities.put("-", 2);
        defaultArities.put("x", 0);
        defaultArities.put("y", 0);
    }

    /**
     * Construct the default symbolic regression language.
     *
     * @param _maxDepth maximum tree depth.
     *
     */
    public Language(int _maxDepth) {

       this(_maxDepth, defaultArities);
    }

    /**
     * Construct a Language.
     *
     * @param _maxDepth maximum tree depth.
     * @param _arities map from Strings representing functions and
     * terminals to their arities.
     */
    public Language(int _maxDepth, HashMap<String, Integer> _arities) {
        maxDepth = _maxDepth;
        arities = _arities;
        terminals = new ArrayList<String>();
        functions = new ArrayList<String>();

        for (String key: arities.keySet()) {
            // System.out.println(key);
            int arity = arities.get(key);
            if (arity == 0) {
                terminals.add(key);
            } else {
                functions.add(key);
            }
            if (arity > maxArity) {
                maxArity = arity;
            }
        }

        Collections.sort(functions);
        Collections.sort(terminals);
        F = functions.size();
        T = terminals.size();

        // PT, the probability of choosing a terminal, defaults to
        // just the proportion of terminals among all node types.
        PT = (float) T / (T + F);

        // "holder" is a special type -- it's the node from which the
        // tree proper hangs. It has arity 1 in that only tree can be
        // attached to it. This is used in some traversal algorithms.
        arities.put("holder", 1);
    }

    public String chooseTerminal(Random rng) {
        return terminals.get(rng.nextInt(terminals.size()));
    }

    public String chooseFunction(Random rng) {
        return functions.get(rng.nextInt(functions.size()));
    }

    // FIXME assumes binary operators only!
    public int countFullTrees() {
        double x = pow(2.0, 4.0);
        return (int) (pow((double) F, pow(2.0, (double) (maxDepth  - 1)))
                      *
                      pow((double) T, pow(2.0, (double) maxDepth)));
    }

    public static void main(String args[]) {
        // If we pass in no arities, it will assume default SR
        // language.
        Language language = new Language(5);
        System.out.println(language.countFullTrees());
    }
}
