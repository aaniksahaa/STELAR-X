package stelarx.gpu;

/**
 * JNI bridge to the CUDA distance-matrix computation kernel.
 *
 * The native implementation uses:
 *   - Euler tour + O(1) sparse-table RMQ for LCA queries
 *   - Δ-tree batching: tree data occupies O(Δ × n log n) GPU VRAM
 *   - B×B pair tiling:  output buffer occupies O(B²) GPU VRAM
 *
 * The full n×n accumulator (dist_sum, cooccur) lives in CPU RAM;
 * the GPU only holds a B×B tile at any moment.
 *
 * Control flow inside the native code:
 *   outer: tree batches  (k/Δ uploads of tree data)
 *   inner: B×B pair tiles (each tile is zeroed, kernel'd, downloaded, CPU-added)
 *
 * GPU VRAM = B² × 12 + Δ × (E_max × 2 + LOG × E_max × 2 + n × 6 + 4)  bytes
 */
public class GPUDistanceMatrix {

    private static volatile boolean loaded    = false;
    private static volatile boolean attempted = false;
    private static volatile String loadError = "not attempted";

    public static synchronized boolean tryLoad() {
        if (attempted) return loaded;
        attempted = true;
        try {
            System.loadLibrary("stelarx_dist");
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
     * Compute the pairwise topological distance matrix on GPU.
     *
     * @param eulerDepths    flat [numTrees × E_max] — Euler tour depths per tree (int16)
     * @param sparseMin      flat [numTrees × LOG × E_max] — sparse table min depths (int16)
     * @param leafDepth      flat [numTrees × n] — depth of each leaf in each tree (int16, -1 absent)
     * @param firstOcc       flat [numTrees × n] — first Euler position of each taxon (-1 absent)
     * @param eulerLen       [numTrees] — actual tour length per tree
     * @param numTrees       k
     * @param n              total taxon count
     * @param E_max          padded Euler tour length (columns in eulerDepths / sparseMin)
     * @param LOG            number of sparse-table levels
     * @param tileSizeB      B — GPU pair tile side length (controls VRAM for output buffer)
     * @param progressInterval  seconds between progress updates (if progressMaxSteps == 0)
     * @param progressMaxSteps  max progress bar prints (0 = use time-interval mode)
     * @param distSumOut     pre-zeroed [n × n] double[] — native code fills this
     * @param cooccurOut     pre-zeroed [n × n] int[]    — native code fills this
     */
    public static native void computeDistancesGPU(
        short[] eulerDepths,
        short[] sparseMin,
        short[] leafDepth,
        int[]   firstOcc,
        int[]   eulerLen,
        int     numTrees,
        int     n,
        int     E_max,
        int     LOG,
        int     tileSizeB,
        double  progressInterval,
        int     progressMaxSteps,
        double[] distSumOut,
        int[]    cooccurOut
    );
}
