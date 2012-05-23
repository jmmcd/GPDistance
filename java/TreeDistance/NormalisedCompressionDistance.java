/*
 * NormalisedCompressionDistance.java
 *
 * Created on August 23 2009
 *
 */
package TreeDistance;

import java.util.zip.*;

/**
 * Calculate the Normalised Compression Distance (NCD), an
 * approximation of the universal (Kolmogorov) similarity metric,
 * between two strings.
 *
 * From "Clustering by Compression", Cilibrasi and Vitanyi
 * http://arXiv.org/abs/cs/0312044v2
 *
 * @author James McDermott
 */
public class NormalisedCompressionDistance {

    public static double ncd(Tree x, Tree y) {
        return ncd(x.toString(), y.toString());
    }

    public static double ncd(String x, String y) {
        // System.out.println("C(x+y): " + C(x+y));
        // System.out.println("C(x): " + C(x));
        // System.out.println("C(y): " + C(y));
        int cx = C(x);
        int cy = C(y);
        int cxy = C(x + y);
        return (cxy - (double) Math.min(cx, cy)) / Math.max(cx, cy);
    }

    // Find the length of a String after compression
    public static int C(String inputString) {
        // Encode a String into bytes
        byte[] input = inputString.getBytes();

        // Compress the bytes
        byte[] output = new byte[input.length + 100];
        Deflater compresser = new Deflater();
        compresser.setInput(input);
        compresser.finish();

        // Return the *length* of the compressed version.
        return compresser.deflate(output);
    }

}
