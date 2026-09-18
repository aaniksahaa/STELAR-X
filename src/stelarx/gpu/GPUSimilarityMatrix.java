package stelarx.gpu;

/**
 * JNI bridge to the CUDA similarity-matrix computation kernel.
 *
 * Implements the validated bridge formula:
 *
 *   same_side_T(a,b)  =  C2(kt − 2)  −  QD_T(a,b)
 *
 *   QD_T(x,y) = ½ · [ (F(x) − F(cx)) + (F(y) − F(cy)) + (cxS−1)·Z + (cyS−1)·Z ]
 *
 * where:
 *   w  = LCA(x, y) in tree T
 *   cx = child of w on the x-side, cy = child of w on the y-side
 *   cxS = s(cx),  cyS = s(cy),   Z = kt − cxS − cyS
 *   F(v) is the path-prefix sum along the root→v path described in
 *   EulerTourBuilder.
 *
 * GPU per-pair query for tree T:
 *   l = min(firstOcc[x], firstOcc[y]),  r = max(...)
 *   k_lvl = floor(log2(r − l + 1)),     l2 = r − 2^k_lvl + 1
 *   read the two compact unsigned-16 argmin positions, compare their depths
 *   with the same left-biased tie rule, then read
 *   (leftChildS, leftChildF, rightChildS, rightChildF) from the base Euler
 *   arrays at the winning LCA INTERMEDIATE position.
 *   If firstOcc[x] ≤ firstOcc[y]:  x's child = left,  y's child = right.
 *   Else: swap.
 *
 * Architecture:
 *   - Δ-tree batching: tree data O(Δ · n · log n) GPU VRAM, capped at
 *     512 MiB by default (configurable via --gpu-sim-vram-cap-mb)
 *   - B×B pair tiling: output tile O(B²) GPU VRAM (B ≈ √(n·k))
 *   - No atomics: thread (da,db) owns pair (a0+da, b0+db) uniquely
 */
public class GPUSimilarityMatrix {

    private static volatile boolean loaded    = false;
    private static volatile boolean attempted = false;
    private static volatile String loadError = "not attempted";

    public static synchronized boolean tryLoad() {
        if (attempted) return loaded;
        attempted = true;
        try {
            System.loadLibrary("stelarx_sim");
            loaded = true;
            loadError = "";
        } catch (UnsatisfiedLinkError | SecurityException e) {
            loaded = false;
            loadError = e.getClass().getSimpleName() + ": " + e.getMessage();
        }
        return loaded;
    }

    public static boolean isLoaded() { return loaded; }
    public static String getLoadError() { return loadError; }

    /**
     * Compute the pairwise similarity matrix on GPU using the bridge formula.
     *
     * @param eulerDepths       flat [numTrees × E_max]            (short) tour depths
     * @param eulerF            flat [numTrees × E_max]            (double) F(node) per pos
     * @param eulerLeftChildS   flat [numTrees × E_max]            (short)  s(leftChild) at intermediates
     * @param eulerLeftChildF   flat [numTrees × E_max]            (double) F(leftChild) at intermediates
     * @param eulerRightChildS  flat [numTrees × E_max]            (short)  s(rightChild) at intermediates
     * @param eulerRightChildF  flat [numTrees × E_max]            (double) F(rightChild) at intermediates
     * @param sparseArgmin      flat [numTrees × LOG × E_max]      (char) unsigned-16
     *                          left-biased argmin Euler position; depth and child
     *                          payloads are fetched from the base Euler arrays
     * @param firstOcc          flat [numTrees × n]                (int)    first tour pos, −1 absent
     * @param eulerLen          [numTrees]                         (int)    actual tour length
     * @param leafCount         [numTrees]                         (int)    kt per tree
     * @param numTrees          k
     * @param n                 total taxon count
     * @param E_max             padded Euler tour length
     * @param LOG               number of sparse-table levels
     * @param tileSizeB         B — GPU pair tile side (0 = auto)
     * @param treeVramCapMiB    maximum tree-batch data allocation in MiB
     * @param progressInterval  seconds between progress updates
     * @param progressMaxSteps  max progress prints (0 = time-interval mode)
     * @param numSumOut         [n × n] accumulator — native adds numerator sums
     * @param denSumOut         [n × n] accumulator — native adds denominator sums
     * @param packedNumOut      segmented packed upper triangle, or null for dense mode
     * @param packedDenOut      segmented packed upper triangle, or null for dense mode
     */
    public static native void computeSimilarityGPU(
        short[]  eulerDepths,
        double[] eulerF,
        short[]  eulerLeftChildS,
        double[] eulerLeftChildF,
        short[]  eulerRightChildS,
        double[] eulerRightChildF,
        char[]   sparseArgmin,
        int[]    firstOcc,
        int[]    eulerLen,
        int[]    leafCount,
        int      numTrees,
        int      n,
        int      E_max,
        int      LOG,
        int      tileSizeB,
        int      treeVramCapMiB,
        double   progressInterval,
        int      progressMaxSteps,
        double[] numSumOut,
        double[] denSumOut,
        double[][] packedNumOut,
        double[][] packedDenOut,
        int packedSegmentShift
    );

    /**
     * Exact large-tour variant. It uses a two-level blocked RMQ, 32-bit Euler
     * depths/child sizes, unsigned-byte in-block argmins, and 32-bit block
     * argmins. The compact method above remains the default for ordinary tours.
     */
    public static native void computeSimilarityGPUWide(
        int[]    eulerDepths,
        double[] eulerF,
        int[]    eulerLeftChildS,
        double[] eulerLeftChildF,
        int[]    eulerRightChildS,
        double[] eulerRightChildF,
        byte[]   microArgmin,
        int[]    macroArgmin,
        int[]    firstOcc,
        int[]    eulerLen,
        int[]    leafCount,
        int      numTrees,
        int      n,
        int      E_max,
        int      microLog,
        int      blockSize,
        int      blockMax,
        int      macroLog,
        int      tileSizeB,
        int      treeVramCapMiB,
        double   progressInterval,
        int      progressMaxSteps,
        double[] numSumOut,
        double[] denSumOut,
        double[][] packedNumOut,
        double[][] packedDenOut,
        int packedSegmentShift
    );
}
