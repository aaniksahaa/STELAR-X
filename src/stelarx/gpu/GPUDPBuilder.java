package stelarx.gpu;

/**
 * JNI bridge to the CUDA cross-tree DP transition search kernel.
 *
 * The native library implements size-binned hash-subtraction search:
 *   For each cluster A in X, for each size sz in [1..|A|/2],
 *   for each cluster B in bin(sz): if hash(A)-hash(B) ∈ X → emit A→B|(A\B).
 *
 * GPU memory is bounded because the output is collected in fixed-size
 * sub-batches (one per sub-batch within each size-sz round).  Sub-batching
 * is transparent to the Java caller; results arrive as one flat int[].
 */
public class GPUDPBuilder {

    private static volatile boolean loaded = false;
    private static volatile boolean loadAttempted = false;
    private static volatile String loadError = "not attempted";

    public static synchronized boolean tryLoad() {
        if (loadAttempted) return loaded;
        loadAttempted = true;
        try {
            System.loadLibrary("stelarx_dp");
            loaded = true;
            loadError = "";
        } catch (UnsatisfiedLinkError | SecurityException e) {
            loaded = false;
            loadError = e.getClass().getSimpleName() + ": " + e.getMessage();
        }
        return loaded;
    }

    public static String getLoadError() { return loadError; }

    /**
     * Find all cross-tree DP transitions on GPU.
     *
     * @param clusterSums  [N*m] cluster sums, clusterSums[c*m+s] = sum of cluster c under seed s
     * @param clusterXors  [N*m] cluster XORs, same layout
     * @param clusterSizes [N]   cluster sizes
     * @param numClusters  N
     * @param m            number of hash seeds
     * @param sortedBySize [N] cluster indices sorted by size ascending
     * @param binStart     [maxSize+2] binStart[sz] = first index in sortedBySize with size >= sz
     * @param maxSize      maximum cluster size present in X
     * @param maxPerRound  max GPU output buffer size per sub-batch (controls VRAM usage)
     * @return int[] = [count, idxA0,idxB0,idxRes0, idxA1,idxB1,idxRes1, ...]
     *         where indices are into the clusterSums/clusterXors/clusterSizes arrays.
     *         Returns null if native library call fails.
     */
    public static native int[] findCrossTreeTransitionsGPU(
        long[]  clusterSums,
        long[]  clusterXors,
        int[]   clusterSizes,
        int     numClusters,
        int     m,
        int[]   sortedBySize,
        int[]   binStart,
        int     maxSize,
        int     maxPerRound,
        double  progressInterval,
        int     progressMaxSteps
    );
}
