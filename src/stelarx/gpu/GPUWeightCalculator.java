package stelarx.gpu;

/**
 * JNI bridge to the CUDA weight-calculation kernel.
 *
 * The native library (libstelarx_weight.so) must be on java.library.path
 * or the directory passed via -Djava.library.path=native/.
 *
 * Call tryLoad() once at startup; it returns false if the .so is missing.
 * After a successful load, computeWeightsGPU() offloads all score computation
 * to the GPU with adaptive split-batching to bound peak VRAM usage.
 *
 * batchSizeHint controls the batching:
 *    0  — auto: the C layer queries cudaMemGetInfo after uploading static data
 *               and picks the largest batch that fits in 75% of free VRAM.
 *   -1  — no batching: all splits in one launch (original single-kernel behaviour).
 *   >0  — manual override: use exactly this many splits per launch.
 */
public class GPUWeightCalculator {

    private static volatile boolean loaded = false;
    private static volatile boolean loadAttempted = false;
    private static volatile String loadError = "not attempted";
    private static volatile Probe cachedProbe;

    /** Result of probing the actual CUDA runtime, independent of nvidia-smi. */
    public record Probe(boolean libraryLoaded, boolean cudaAvailable,
                        String deviceName, int deviceCount,
                        int computeMajor, int computeMinor,
                        int driverVersion, int runtimeVersion,
                        long freeMiB, long totalMiB, String detail) {}

    /** Try to load the native library; returns true on success. */
    public static synchronized boolean tryLoad() {
        if (loadAttempted) return loaded;
        loadAttempted = true;
        try {
            // Distinct name prevents an older STELAR-X quartet kernel from being
            // loaded accidentally into the rooted STELAR-X scoring pipeline.
            System.loadLibrary("stelarx_weight");
            loaded = true;
            loadError = "";
        } catch (UnsatisfiedLinkError | SecurityException e) {
            // Not a fatal error — caller falls back to CPU path
            loaded = false;
            loadError = e.getClass().getSimpleName() + ": " + e.getMessage();
        }
        return loaded;
    }

    public static boolean isLoaded() { return loaded; }
    public static String getLoadError() { return loadError; }

    /**
     * Load the backend and ask CUDA itself whether a usable device/driver exists.
     * The result is cached because probing initializes the CUDA runtime.
     */
    public static synchronized Probe probe() {
        if (cachedProbe != null) return cachedProbe;
        if (!tryLoad()) {
            cachedProbe = new Probe(false, false, "", 0, 0, 0, 0, 0, 0, 0,
                "native CUDA backend could not be loaded: " + loadError);
            return cachedProbe;
        }
        try {
            String status = queryGPUStatus();
            cachedProbe = parseProbe(status);
        } catch (Throwable t) {
            cachedProbe = new Probe(true, false, "", 0, 0, 0, 0, 0, 0, 0,
                "CUDA probe failed: " + t.getClass().getSimpleName() + ": " + t.getMessage());
        }
        return cachedProbe;
    }

    private static Probe parseProbe(String status) {
        if (status == null || status.isBlank()) {
            return new Probe(true, false, "", 0, 0, 0, 0, 0, 0, 0,
                "CUDA probe returned no status");
        }
        java.util.Map<String,String> kv = new java.util.HashMap<>();
        String[] fields = status.split(";", -1);
        for (int i = 1; i < fields.length; i++) {
            int eq = fields[i].indexOf('=');
            if (eq > 0) kv.put(fields[i].substring(0, eq), fields[i].substring(eq + 1));
        }
        boolean ok = fields[0].equals("OK");
        return new Probe(true, ok,
            kv.getOrDefault("name", ""), intValue(kv, "devices"),
            intValue(kv, "ccMajor"), intValue(kv, "ccMinor"),
            intValue(kv, "driver"), intValue(kv, "runtime"),
            longValue(kv, "freeMiB"), longValue(kv, "totalMiB"),
            kv.getOrDefault("detail", status));
    }

    private static int intValue(java.util.Map<String,String> kv, String key) {
        try { return Integer.parseInt(kv.getOrDefault(key, "0")); }
        catch (NumberFormatException e) { return 0; }
    }

    private static long longValue(java.util.Map<String,String> kv, String key) {
        try { return Long.parseLong(kv.getOrDefault(key, "0")); }
        catch (NumberFormatException e) { return 0L; }
    }

    /**
     * Compute 2*score for every split on the GPU (prefix-sum tree-DP) with
     * adaptive split-batching.
     *
     * One thread block per split; the block loops over every gene tree, builds
     * per-tree leaf prefix sums for both sides of the split in shared memory,
     * then derives the child intersections for each rooted internal partition
     * in O(1). Binary and polytomous nodes use the same rooted-triplet objective.
     *
     * Static data (orderings, invIndex, node CSR) is uploaded to the GPU once.
     * Splits are streamed in adaptive batches; scores are accumulated into the
     * host result array.  This bounds peak VRAM at:
     *
     *   O(numGpuTrees × numTaxa)   [orderings + invIndex, permanent]
     * + O(totalNodes)              [node CSR, permanent]
     * + batchSize × 48 B           [current split batch + score slice]
     *
     * Per-block transient working set is 2·(maxLeafCount+1) ints of shared memory.
     *
     * @param splits         flat int array, numSplits × 10
     *                       [aTree,aLo,aHi,aComp,aSize, bTree,bLo,bHi,bComp,bSize]
     * @param nodeData       flat int array, numUnique × 3  [lo, mid, hi] (exemplar interval)
     * @param nodeFreq       flat int array, numUnique  (frequency of each unique binary partition)
     * @param nodeOffset     flat int array, numPartTrees + 1  (CSR row pointers, bucket by exemplar)
     * @param partLeafCount  flat int array, numPartTrees  (leaf count L per gene tree)
     * @param orderings      flat int array, numGpuTrees × numTaxa
     * @param invIndex       flat int array, numGpuTrees × numTaxa
     * @param numSplits      number of candidate splits
     * @param numPartTrees   number of gene trees contributing rooted child partitions
     * @param partTreeOffset orderings/invIndex slot offset for gene trees
     *                       (0, or numClusterTrees when autocomplete is active)
     * @param maxLeafCount   max leaf count over the gene trees (shared-mem sizing)
     * @param numGpuTrees    total orderings/invIndex slots (for VRAM accounting)
     * @param numTaxa        total taxon count (registry size)
     * @param batchSizeHint  0=auto, -1=no batching, >0=exact batch size
     * @param vramFraction   fraction of free VRAM to use when batchSizeHint==0
     * @param scoreMode      accumulator/transport selector:
     *                       0 = LONG   — exact 64-bit integer 2·score per slot;
     *                       1 = DOUBLE — IEEE-754 bit pattern of the 2·score per
     *                                    slot (decode with {@link Double#longBitsToDouble});
     *                       2 = INT128 — exact 128-bit 2·score as TWO longs per
     *                                    split: result[2*i]=low (unsigned),
     *                                    result[2*i+1]=high (signed).
     * @return for LONG/DOUBLE: long[numSplits]; for INT128: long[2*numSplits];
     *         or null if the GPU path is infeasible (caller falls back to CPU)
     */
    public static native long[] computeWeightsGPU(
        int[] splits,
        int[] splitRangeMeta,   // numSplits*4: [aRngOff,aRngCnt,bRngOff,bRngCnt]; cnt 0 = single-range
        int[] rangeData,        // flat [lo,hi] pairs for multi-range split sides (resident)
        int[] nodeData,
        int[] nodeFreq,
        int[] nodeOffset,
        int[] partLeafCount,
        // Polytomy (d>3) CSR — empty (polyTreeOffset all-zero, others length 0/1)
        // when there are no polytomous partitions ⇒ kernel poly loop is a no-op.
        int[] polyTreeOffset,   // numPartTrees+1: CSR row pointers (poly nodes bucketed by tree)
        int[] polyBoundOffset,  // numPoly+1: range into polyBounds (length d) per poly node
        int[] polyBounds,       // concatenated boundary lists b[0..d-1] (child i = [b[i],b[i+1]))
        int[] polyFreq,         // numPoly: occurrence count per unique poly partition
        int[] orderings,
        int[] invIndex,
        int numSplits,
        int numPartTrees,
        int partTreeOffset,
        int maxLeafCount,
        int numGpuTrees,
        int numTaxa,
        int batchSizeHint,
        double vramFraction,
        int scoreMode,
        double progressIntervalSec     // -1 = auto; else override the progress cadence
    );

    /**
     * Legacy "smaller-side traversal" weight calculation (no prefix sums).
     *
     * One CUDA thread per split, zero per-thread working state: each thread loops
     * over every deduplicated rooted child partition and computes the required
     * intersections by walking the smaller of the two ranges element-by-element
     * (looking taxa up in invIndex).  No prefix-sum arrays are built or allocated,
     * so device memory is just the static parts/orderings/invIndex plus the split
     * batch — there is no O(L) prefix working set.
     *
     * Static data (parts, orderings, invIndex) is uploaded once; splits stream in
     * adaptive batches exactly like the prefix-sum path.
     *
     * @param splits     flat int array, numSplits × 10
     *                   [aTree,aLo,aHi,aComp,aSize, bTree,bLo,bHi,bComp,bSize]
     * @param parts      flat int array, numParts × 8 (deduplicated binary partitions)
     *                   [treeIdx, lo1, hi1, lo2, hi2, sz1, sz2, frequency]
     * @param orderings  flat int array, numGpuTrees × numTaxa
     * @param invIndex   flat int array, numGpuTrees × numTaxa
     * @param numSplits  number of candidate splits
     * @param numParts   number of unique binary rooted child partitions
     * @param numGpuTrees total orderings/invIndex slots
     * @param numTaxa    total taxon count (registry size)
     * @param totalN     total taxon count (same as numTaxa, used for sizeC)
     * @param batchSizeHint 0=auto, -1=no batching, >0=exact batch size
     * @param vramFraction  fraction of free VRAM to use when batchSizeHint==0
     * @param scoreMode     accumulator/transport selector (see computeWeightsGPU):
     *                      0=LONG, 1=DOUBLE (bit pattern), 2=INT128 (low,high pair).
     * @return for LONG/DOUBLE: long[numSplits]; for INT128: long[2*numSplits];
     *         or null on failure
     */
    public static native long[] computeWeightsSmallerSideGPU(
        int[] splits,
        int[] splitRangeMeta,   // numSplits*4: [aRngOff,aRngCnt,bRngOff,bRngCnt]; cnt 0 = single-range
        int[] rangeData,        // flat [lo,hi] pairs for multi-range split sides (resident)
        int[] parts,
        // Polytomy (d>3) CSR — empty when no polytomous partitions.
        int[] ssPolyMeta,       // numPolyParts*2: [treeIdx(+offset), freq]
        int[] ssPolyBoundOffset,// numPolyParts+1: range into ssPolyBounds (length d) per poly node
        int[] ssPolyBounds,     // concatenated boundary lists b[0..d-1]
        int[] orderings,
        int[] invIndex,
        int numSplits,
        int numParts,
        int numPolyParts,
        int numGpuTrees,
        int numTaxa,
        int totalN,
        int batchSizeHint,
        double vramFraction,
        int scoreMode,
        double progressIntervalSec     // -1 = auto; else override the progress cadence
    );

    /**
     * Bitset weight calculation (low-taxa fast path).
     *
     * Every cluster and every gene-tree part is materialized on the host as a
     * global-taxon bitset of W = ceil(numTaxa/64) 64-bit words.  One CUDA thread
     * per split loads its A/B cluster bitsets (indexed by a per-split cluster id
     * into the resident cluster pool) and, for each unique part, computes each core
     * intersection as popcount(A & M) over W words — no orderings/invIndex or prefix
     * arrays are used in the kernel. The child-only QI formula is identical
     * to the other two methods, so scores are bit-identical.
     *
     * Resident (uploaded once): clusterBits, partM1/partM2 + frequencies,
     * poly CSR.  Per-batch (streamed like the other paths): splits (4 ints/split) +
     * scores.  The batching / VRAM logic matches computeWeightsSmallerSideGPU.
     *
     * @param splits        curBatch source: numSplits × 4  [aCid, bCid, aSize, bSize]
     * @param clusterBits   numClusters × W  longs (global-taxon bitsets; cid 0 = empty)
     * @param partM1        numParts × W  longs (binary part M1 bitset)
     * @param partM2        numParts × W  longs (binary part M2 bitset)
     * @param partFreq      numParts occurrence counts
     * @param polyMeta      numPoly × 2 ints [childCount, freq]
     * @param polyChildOffset numPoly + 1 CSR row pointers into polyChildBits
     * @param polyChildBits ΣchildCount × W longs
     * @param numSplits     number of candidate splits
     * @param numClusters   number of distinct cluster bitsets in the pool
     * @param numParts      number of unique binary (d==3) parts
     * @param numPoly       number of unique polytomous (d>3) parts
     * @param wordsPerSet   W = ceil(numTaxa/64)
     * @param numTaxa       total taxon count (= totalN for sizeC)
     * @param batchSizeHint 0=auto, -1=no batching, >0=exact batch size
     * @param vramFraction  fraction of free VRAM to use when batchSizeHint==0
     * @param scoreMode     0=LONG, 1=DOUBLE (bit pattern), 2=INT128 (low,high pair)
     * @return for LONG/DOUBLE: long[numSplits]; for INT128: long[2*numSplits]; or null on failure
     */
    public static native long[] computeWeightsBitsetGPU(
        int[]  splits,
        long[] clusterBits,
        long[] partM1,
        long[] partM2,
        int[]  partFreq,
        int[]  polyMeta,
        int[]  polyChildOffset,
        long[] polyChildBits,
        int numSplits,
        int numClusters,
        int numParts,
        int numPoly,
        int wordsPerSet,
        int numTaxa,
        int batchSizeHint,
        double vramFraction,
        int scoreMode,
        double progressIntervalSec
    );

    /**
     * Simple-tree-walk weight calculation (many-candidate fast path).
     *
     * One CUDA thread per split walks a resident flat postorder token stream of all
     * gene trees sequentially, maintaining a small per-thread stack of
     * (|node∩A|,|node∩B|) pairs. Every rooted internal node,
     * including the gene-tree root, is scored from its actual children using the
     * rooted-triplet objective (bit-identical to the other methods). No prefix
     * arrays, no dedup, no cross-tree
     * parallelism — a lean kernel that wins when the candidate set is huge.
     *
     * Resident (uploaded once): clusterBits (A/B pool), nodeStream and
     * treeNodeOffset. Per-batch (streamed): splits (4 ints/split) + scores.
     *
     * @param splits         numSplits × 4  [aCid, bCid, aSize, bSize]
     * @param clusterBits    numClusters × W  longs (global-taxon bitsets; cid 0 = empty)
     * @param nodeStream     flat postorder tokens: leaf = taxon id (≥0), internal = -childCount
     * @param treeNodeOffset numTrees + 1  CSR row pointers into nodeStream
     * @param numSplits      number of candidate splits
     * @param numClusters    number of distinct cluster bitsets in the pool
     * @param numTrees       number of gene trees
     * @param wordsPerSet    W = ceil(numTaxa/64)
     * @param numTaxa        total taxon count (= totalN for sizeC)
     * @param maxFrontier    measured maximum postorder stack entries over all trees
     * @param batchSizeHint  0=auto, -1=no batching, >0=exact batch size
     * @param vramFraction   fraction of free VRAM to use when batchSizeHint==0
     * @param scoreMode      0=LONG, 1=DOUBLE (bit pattern), 2=INT128 (low,high pair)
     * @return for LONG/DOUBLE: long[numSplits]; for INT128: long[2*numSplits];
     *         or null if infeasible (e.g. maxFrontier exceeds the compiled GPU
     *         stack cap) — caller falls back to CPU
     */
    public static native long[] computeWeightsTreeWalkGPU(
        int[]  splits,
        long[] clusterBits,
        int[]  nodeStream,
        int[]  treeNodeOffset,
        int numSplits,
        int numClusters,
        int numTrees,
        int wordsPerSet,
        int numTaxa,
        int maxFrontier,
        int batchSizeHint,
        double vramFraction,
        int scoreMode,
        double progressIntervalSec
    );

    /**
     * Query GPU free and total VRAM via cudaMemGetInfo.
     * Returns long[2] = {freeMiB, totalMiB}, or null if unavailable.
     */
    public static native long[] queryVRAMMiB();

    /** Machine-readable CUDA/device probe used for auto-selection and diagnostics. */
    private static native String queryGPUStatus();
}
