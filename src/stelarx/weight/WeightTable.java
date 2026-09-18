package stelarx.weight;

import stelarx.Config;
import stelarx.Logging;
import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.dp.BipartitionSplit;
import stelarx.dp.DPTable;
import stelarx.gpu.GPUWeightCalculator;
import stelarx.partition.Partition;
import stelarx.partition.PartitionTable;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.Int128;
import stelarx.util.ProgressBar;
import stelarx.util.Threading;

import java.util.*;

/**
 * Precomputed STELAR-X rooted-triplet scores for every candidate child split.
 *
 * For a candidate rooted child split (A | B), and a rooted gene-tree child
 * partition (M1 | M2), compute:
 *
 *   a0=|A∩M1|, a1=|A∩M2|, b0=|B∩M1|, b1=|B∩M2|
 *
 *   2w = a0*b1*(a0+b1-2) + a1*b0*(a1+b0-2)
 *
 *   score(A|B) = sum over unique rooted gene-tree partitions P:
 *                P.frequency * w
 *
 * The doubled expression is always even, so each weight is an exact
 * non-negative integer and
 * so scores are stored as non-negative longs.
 *
 * Execution path selection:
 *   GPU  — when --gpu is set and libstelarx_weight.so is loadable.
 *   CPU  — otherwise (multi-threaded via Threading.processRangeParallel).
 */
public class WeightTable {

    /**
     * Numeric type used for scores. Decided once per run from the problem size:
     * LONG below the overflow threshold; above it, INT128 (exact, default) or
     * DOUBLE (approximate) per {@link Config#getLargeScoreType()}.
     */
    public enum Mode { LONG, DOUBLE, INT128 }

    // DPTable already owns every BipartitionSplit. Assign each scored split a
    // dense index and keep only the mode-matching value array here. This avoids
    // a duplicate HashMap plus boxed Long/Double values at peak residency.
    private long[]   scores  = new long[0];
    private double[] scoresD = new double[0];
    private Int128[] scoresI = new Int128[0];
    private BipartitionSplit[] scoredSplits = new BipartitionSplit[0];
    private final int n;   // total taxa

    private final Mode mode;
    private final boolean useDouble;   // mode == DOUBLE
    private final boolean useInt128;   // mode == INT128

    // stats (LONG path)
    private long maxScore;
    private long totalScore;
    // stats (DOUBLE path)
    private double maxScoreD = Double.NEGATIVE_INFINITY;
    private double totalScoreD;
    // stats (INT128 path)
    private Int128 maxScoreI   = null;          // null = -inf sentinel
    private Int128 totalScoreI = Int128.ZERO;

    // -------------------------------------------------------------------------

    /**
     * @param clusterTrees  completed gene trees — used for cluster exemplar position lookups
     *                      (Cluster.treeIndex refers into this list)
     * @param partTrees     original (pre-completion) gene trees — used for rooted-triplet
     *                      scoring (Partition.treeIndex refers into this list).
     *                      When --autocomplete-incomplete-gene-trees is NOT active,
     *                      partTrees == clusterTrees (same reference) and behavior is unchanged.
     */
    public WeightTable(DPTable dpTable, PartitionTable partTable,
                       ClusterTable clusterTable,
                       List<Tree> clusterTrees, List<Tree> partTrees) {
        long t0 = System.nanoTime();
        this.n = clusterTable.getAllTaxaHash().size;

        // ── Numeric-precision decision (LONG vs DOUBLE/INT128 accumulation) ───
        // The exact triplet score grows as O(genes · n^3) and overflows a signed
        // 64-bit integer for very large taxon sets.  Below the overflow threshold
        // we keep exact LONG.  Above it we use INT128 (exact, default) or DOUBLE
        // (approximate) per the configured large-score type.
        int numGenes = partTrees.size();
        if (needsDoubleAccumulation(n, numGenes)) {
            this.mode = (Config.getInstance().getLargeScoreType() == Config.LargeScoreType.DOUBLE)
                        ? Mode.DOUBLE : Mode.INT128;
        } else {
            this.mode = Mode.LONG;
        }
        this.useDouble = (mode == Mode.DOUBLE);
        this.useInt128 = (mode == Mode.INT128);
        logAccumulationDecision(n, numGenes, mode);

        // Collect all unique splits from DPTable into an indexed list.  Optionally
        // prune splits whose parent cluster is unreachable from the DP root — the
        // top-down inference DP never scores those, so skipping them is
        // result-preserving (see DPTable.reachableClusters).
        List<BipartitionSplit> splitList = new ArrayList<>();
        if (Config.getInstance().isPruneUnreachableSplits()) {
            Set<ClusterHash> reachable = dpTable.reachableClusters();
            int totalClusters = 0, keptClusters = 0, totalSplits = 0;
            for (var entry : dpTable.entries()) {
                totalClusters++;
                int sz = entry.getValue().size();
                totalSplits += sz;
                if (reachable.contains(entry.getKey())) {
                    splitList.addAll(entry.getValue());
                    keptClusters++;
                }
            }
            int kept = splitList.size();
            Logging.info("DP reachability prune: clusters %d/%d reachable (%.1f%%), "
                + "splits %d/%d need to be scored (%.1f%%) — %d unreachable splits skipped",
                keptClusters, totalClusters, 100.0 * keptClusters / Math.max(1, totalClusters),
                kept, totalSplits, 100.0 * kept / Math.max(1, totalSplits), totalSplits - kept);
        } else {
            for (var entry : dpTable.entries()) splitList.addAll(entry.getValue());
        }
        int numSplits = splitList.size();

        scoredSplits = splitList.toArray(BipartitionSplit[]::new);
        for (int i = 0; i < numSplits; i++) scoredSplits[i].assignScoreIndex(i);

        // Nothing below is meaningful for an empty launch, and the native auto-
        // batching code necessarily divides by the resolved batch size.  This can
        // occur for a malformed or externally constructed DP table; return a valid
        // empty weight table instead of entering JNI with numSplits == 0.
        if (numSplits == 0) {
            Logging.info("Weight table: no reachable splits to score");
            return;
        }

        // Per-split score buffers — exactly one is non-null, matching the mode.
        long[]   scoreArray  = (mode == Mode.LONG)   ? new long[numSplits]    : null;  // exact integer
        double[] scoreArrayD = (mode == Mode.DOUBLE) ? new double[numSplits]  : null;  // floating point
        Int128[] scoreArrayI = (mode == Mode.INT128) ? new Int128[numSplits]  : null;  // exact 128-bit

        // When clusterTrees != partTrees, the GPU path packs both sets of
        // orderings/invIndex into a combined array (slots 0..C-1 = cluster trees,
        // slots C..C+P-1 = partition trees, C = clusterTrees.size()) and offsets
        // partition tree indices by C.  numGpuTrees is the combined count.
        //   - autocomplete inference: C == P == k  (completed vs original gene trees)
        //   - score-only mode:        C == 1, P == numGeneTrees  (single species tree)
        // The general C + P form is byte-identical to the old 2k for the C == P case.
        boolean splitTrees = (clusterTrees != partTrees);
        int numGpuTrees = splitTrees ? clusterTrees.size() + partTrees.size() : clusterTrees.size();

        boolean useGPU = (Config.getInstance().getComputeMode() == Config.ComputeMode.GPU)
                         && GPUWeightCalculator.tryLoad();

        // Rooted polytomies are handled natively by every CUDA intersection path.
        // All representations contain actual children only; no outside-node set
        // is hashed, packed, or intersected.

        // Multi-range clusters (consensus emission bridge) are handled fully ON GPU
        // via the two-tier range-CSR: buildSplitRangeData packs each multi-range split
        // side's ranges (single-range sides carry count 0 → the byte-identical fast
        // path), and both kernels sum the intersection over a side's ranges
        // No CPU correction needed.

        Config.WeightIntersectionMethod method = Config.getInstance().getWeightIntersectionMethod();
        boolean bitset   = (method == Config.WeightIntersectionMethod.BITSET);
        boolean treeWalk = (method == Config.WeightIntersectionMethod.SIMPLE_TREE_WALK);

        if (useGPU && treeWalk) {
            boolean ok = computeScoresGPUTreeWalk(splitList, clusterTable,
                                                  clusterTrees, partTrees,
                                                  scoreArray, scoreArrayD, scoreArrayI);
            if (!ok) {
                Logging.info("GPU tree-walk weight path infeasible, falling back to CPU tree-walk");
                computeScoresCPUTreeWalk(splitList, clusterTable, clusterTrees, partTrees,
                                         scoreArray, scoreArrayD, scoreArrayI);
            }
        } else if (useGPU && bitset) {
            boolean ok = computeScoresGPUBitset(splitList, partTable, clusterTable,
                                                clusterTrees, partTrees,
                                                scoreArray, scoreArrayD, scoreArrayI);
            if (!ok) {
                Logging.info("GPU bitset weight path infeasible, falling back to CPU bitset");
                computeScoresCPUBitset(splitList, partTable, clusterTable,
                                       clusterTrees, partTrees,
                                       scoreArray, scoreArrayD, scoreArrayI);
            }
        } else if (useGPU) {
            Config cfg = Config.getInstance();
            boolean smallerSide = (cfg.getWeightIntersectionMethod()
                                   == Config.WeightIntersectionMethod.SMALLER_SIDE_TRAVERSAL);

            // Each path keeps its own resident-data representation.  PREFIX_SUM builds
            // the deduplicated node CSR (and the per-block prefix working memory);
            // SMALLER_SIDE_TRAVERSAL builds none of that — it streams the parts and
            // walks the smaller side per intersection, with zero per-thread state.
            NodeCSR csr = smallerSide ? null : buildDedupNodeCSR(partTable, partTrees);

            // Resident data memory (for the vram-control-factor sizing only).
            long orderingMem = (long) numGpuTrees * n * 2 * Integer.BYTES; // orderings + invIndex
            long modeDataMem; String modeDataDesc;
            if (smallerSide) {
                modeDataMem  = (long) partTable.size() * 8 * Integer.BYTES; // parts
                modeDataDesc = "parts";
            } else {
                modeDataMem  = (long) csr.nodeData.length      * Integer.BYTES
                             + (long) csr.nodeFreq.length      * Integer.BYTES
                             + (long) csr.nodeOffset.length    * Integer.BYTES
                             + (long) csr.partLeafCount.length * Integer.BYTES;
                modeDataDesc = "nodeCSR";
            }

            // Resolve batchSizeHint
            //   Priority: no-batch  >  gpu-batches  >  gpu-batch-size
            //           > gpu-vram-control-factor (explicit)  >  auto (gpu-vram-occupancy-factor)
            //   -1  = no batching (single launch)
            //    0  = auto: native queries free VRAM and computes batch size itself
            //   >0  = exact splits-per-batch resolved here; native uses it directly
            int batchSizeHint;
            String batchDesc;
            if (!cfg.isGpuBatch()) {
                batchSizeHint = -1;
                batchDesc = "off (single launch)";
            } else if (cfg.getGpuNumBatches() > 0) {
                int N = cfg.getGpuNumBatches();
                batchSizeHint = (numSplits + N - 1) / N;
                batchDesc = N + " batches → batchSize=" + batchSizeHint;
            } else if (cfg.getGpuBatchSize() > 0) {
                batchSizeHint = cfg.getGpuBatchSize();
                batchDesc = "explicit batchSize=" + batchSizeHint;
            } else if (cfg.isGpuVramControlFactorSet()) {
                // Manual resident-relative sizing:  mem(batch) = F × mem(resident)
                double F           = cfg.getGpuVramControlFactor();
                long   residentMem = modeDataMem + orderingMem;
                long   batchMem    = (long)(F * residentMem);
                long   perSplit    = 10L * Integer.BYTES + Long.BYTES;              // 48 B/split
                batchSizeHint      = (int) Math.max(1, Math.min(numSplits, batchMem / perSplit));
                int numBatches     = (numSplits + batchSizeHint - 1) / batchSizeHint;
                batchDesc = String.format(
                    "vram-control-factor=%.3f  resident=%.1f MB (%s=%.1f orderings=%.1f)  batch=%.1f MB  → %d batches",
                    F, residentMem / 1e6, modeDataDesc, modeDataMem / 1e6, orderingMem / 1e6, batchMem / 1e6, numBatches);
            } else {
                // Default: auto — pass 0 to native; native queries free VRAM after static upload
                // and computes batchSize = floor(freeVRAM * vramFraction / 48 B)
                batchSizeHint = 0;
                batchDesc = String.format("auto (free-VRAM adaptive, occupancy=%.0f%%)",
                    cfg.getGpuVramFraction() * 100);
            }

            boolean ok;
            if (smallerSide) {
                Logging.info("Weight table: GPU path (smaller-side traversal)  splits=%d  uniqueParts=%d  trees=%d  batching=%s",
                    numSplits, partTable.size(), partTrees.size(), batchDesc);
                ok = computeScoresGPUSmallerSide(splitList, partTable, clusterTable,
                                                 clusterTrees, partTrees, numGpuTrees,
                                                 scoreArray, scoreArrayD, scoreArrayI,
                                                 batchSizeHint, cfg.getGpuVramFraction());
            } else {
                Logging.info("Weight table: GPU path (prefix-sum tree-DP)  splits=%d  uniqueParts=%d  trees=%d  maxLeaf=%d  batching=%s",
                    numSplits, csr.totalNodes, partTrees.size(), csr.maxLeafCount, batchDesc);
                ok = computeScoresGPUPrefixSum(splitList, csr, clusterTable, clusterTrees, partTrees,
                                               numGpuTrees, scoreArray, scoreArrayD, scoreArrayI,
                                               batchSizeHint, cfg.getGpuVramFraction());
            }
            if (!ok) {
                Logging.info("GPU weight path infeasible (e.g. shared-memory limit), falling back to CPU");
                computeScoresCPU(splitList, partTable.entries(), clusterTable,
                                 clusterTrees, partTrees, scoreArray, scoreArrayD, scoreArrayI);
            }
        } else {
            if (Config.getInstance().getComputeMode() == Config.ComputeMode.GPU) {
                Logging.info("GPU library not available, falling back to CPU");
            }
            if (treeWalk) {
                computeScoresCPUTreeWalk(splitList, clusterTable, clusterTrees, partTrees,
                                         scoreArray, scoreArrayD, scoreArrayI);
            } else if (bitset) {
                computeScoresCPUBitset(splitList, partTable, clusterTable,
                                       clusterTrees, partTrees, scoreArray, scoreArrayD, scoreArrayI);
            } else {
                computeScoresCPU(splitList, partTable.entries(), clusterTable,
                                 clusterTrees, partTrees, scoreArray, scoreArrayD, scoreArrayI);
            }
        }

        long ms;
        if (mode == Mode.INT128) {
            scoresI = scoreArrayI;
            for (int i = 0; i < numSplits; i++) {
                Int128 s = scoreArrayI[i];
                if (maxScoreI == null || s.compareTo(maxScoreI) > 0) maxScoreI = s;
                totalScoreI = totalScoreI.add(s);
            }
            ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("Weight table: %d splits scored [INT128], maxScore=%s, totalScore=%s in %d ms",
                numSplits, maxScoreI, totalScoreI, ms);
        } else if (mode == Mode.DOUBLE) {
            scoresD = scoreArrayD;
            for (int i = 0; i < numSplits; i++) {
                double s = scoreArrayD[i];
                if (s > maxScoreD) maxScoreD = s;
                totalScoreD += s;
            }
            ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("Weight table: %d splits scored [DOUBLE], maxScore=%.6e, totalScore=%.6e in %d ms",
                numSplits, maxScoreD, totalScoreD, ms);
        } else {
            scores = scoreArray;
            for (int i = 0; i < numSplits; i++) {
                if (scoreArray[i] > maxScore) maxScore = scoreArray[i];
                totalScore += scoreArray[i];
            }
            ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("Weight table: %d splits scored [LONG], maxScore=%d, totalScore=%d in %d ms",
                numSplits, maxScore, totalScore, ms);
        }
    }

    // -------------------------------------------------------------------------
    // Numeric-precision decision: exact LONG vs floating-point DOUBLE
    // -------------------------------------------------------------------------

    /**
     * Decide whether weight scores must be accumulated as {@code double} to avoid
     * 64-bit integer overflow.
     *
     * <p>For a candidate split scored against {@code numGenes} gene trees, the
     * exact doubled rooted-triplet score is bounded by
     * {@code numGenes · 2 · C(n,3) ≈ numGenes · n^3 / 3}. We compare this
     * estimate against {@code Long.MAX_VALUE / 8} (an 8× safety margin that also
     * covers intermediate {@code frequency · doubled-weight} products and partial sums). When
     * the estimate exceeds that bound, exact {@code long} arithmetic would wrap
     * around (producing the notorious negative scores), so we switch to
     * {@code double}.
     *
     * <p>Overridable for testing via environment variables
     * {@code STELARX_WEIGHT_FORCE_DOUBLE} / {@code STELARX_WEIGHT_FORCE_LONG}.
     */
    static boolean needsDoubleAccumulation(int n, int numGenes) {
        if (System.getenv("STELARX_WEIGHT_FORCE_DOUBLE") != null) return true;
        if (System.getenv("STELARX_WEIGHT_FORCE_LONG")   != null) return false;
        return estimatedMaxTwoScore(n, numGenes) > longSafeBound();
    }

    /** Maximum possible doubled rooted-triplet score: numGenes · 2 · C(n,3). */
    private static double estimatedMaxTwoScore(int n, int numGenes) {
        double nn = (double) n;
        return (double) numGenes * nn * (nn - 1.0) * (nn - 2.0) / 3.0;
    }

    /** Long-safe bound with an 8× margin for intermediate products/sums. */
    private static double longSafeBound() {
        return (double) Long.MAX_VALUE / 8.0;   // ≈ 1.153e18
    }

    /** scoreMode int passed to the native kernel: 0=LONG, 1=DOUBLE, 2=INT128. */
    private int nativeScoreMode() {
        switch (mode) {
            case DOUBLE: return 1;
            case INT128: return 2;
            default:     return 0;
        }
    }

    /** Emit a prominent, human-readable log line stating the chosen score type and why. */
    private static void logAccumulationDecision(int n, int numGenes, Mode mode) {
        double est  = estimatedMaxTwoScore(n, numGenes);
        double safe = longSafeBound();
        switch (mode) {
            case INT128 -> Logging.info(
                "Weight accumulation: INT128 (exact 128-bit integer)  "
                + "[taxa=%d, genes=%d, est. max 2·score ≈ %.2e exceeds long-safe %.2e]  "
                + "— switched to avoid 64-bit integer overflow; scores remain exact "
                + "(full-rate integer math, no FP64 penalty).  Override with "
                + "--large-n-score-type double.", n, numGenes, est, safe);
            case DOUBLE -> Logging.info(
                "Weight accumulation: DOUBLE (64-bit floating point, ~15-16 significant digits)  "
                + "[taxa=%d, genes=%d, est. max 2·score ≈ %.2e exceeds long-safe %.2e]  "
                + "— switched to avoid 64-bit integer overflow; scores are approximate but "
                + "topologically equivalent.", n, numGenes, est, safe);
            default -> Logging.info(
                "Weight accumulation: LONG (exact 64-bit integer)  "
                + "[taxa=%d, genes=%d, est. max 2·score ≈ %.2e within long-safe %.2e].",
                n, numGenes, est, safe);
        }
    }

    // -------------------------------------------------------------------------
    // GPU path
    // -------------------------------------------------------------------------

    /**
     * Flatten all data to primitive arrays, call the CUDA kernel via JNI,
     * and write results (score = twoScore/2) into scoreArray.
     */
    /**
     * GPU weight calculation.
     *
     * When clusterTrees == partTrees (no autocomplete), the orderings/invIndex array has
     * numGpuTrees = k slots (indices 0..k-1) and partition treeIndex values are unchanged.
     *
     * When clusterTrees != partTrees (autocomplete active), numGpuTrees = 2k:
     *   slots 0..k-1   → completed tree orderings/invIndex  (for cluster lookups)
     *   slots k..2k-1  → original tree orderings/invIndex   (for gene-tree/partition lookups)
     * Partition treeIndex values are stored as (p.treeIndex + k) so the kernel naturally
     * reads from the original-tree half of the combined array — no kernel changes needed.
     */
    private boolean computeScoresGPUPrefixSum(List<BipartitionSplit> splitList,
                                               NodeCSR csr,
                                               ClusterTable clusterTable,
                                               List<Tree> clusterTrees,
                                               List<Tree> partTrees,
                                               int numGpuTrees,
                                               long[] scoreArray,
                                               double[] scoreArrayD,
                                               Int128[] scoreArrayI,
                                               int batchSizeHint,
                                               double vramFraction) {
        int numSplits      = splitList.size();
        int numPartTrees   = partTrees.size();
        boolean splitTrees = (clusterTrees != partTrees);
        int partTreeOffset = splitTrees ? clusterTrees.size() : 0;

        int[] splitsData = buildSplitsData(splitList, clusterTable);
        int[][] rng      = buildSplitRangeData(splitList, clusterTable);
        int[] splitRangeMeta = rng[0], rangeData = rng[1];
        int[][] oi       = buildOrderingsInvIndex(clusterTrees, partTrees, numGpuTrees);
        int[] orderings  = oi[0], invIndex = oi[1];

        long t1 = System.nanoTime();
        long[] twoScores = GPUWeightCalculator.computeWeightsGPU(
            splitsData, splitRangeMeta, rangeData,
            csr.nodeData, csr.nodeFreq, csr.nodeOffset, csr.partLeafCount,
            csr.polyTreeOffset, csr.polyBoundOffset, csr.polyBounds, csr.polyFreq,
            orderings, invIndex,
            numSplits, numPartTrees, partTreeOffset, csr.maxLeafCount,
            numGpuTrees, n,
            batchSizeHint, vramFraction, nativeScoreMode(),
            Config.getInstance().getGpuProgressIntervalSec());
        long gpuMs = (System.nanoTime() - t1) / 1_000_000;

        splitsData = null; orderings = null; invIndex = null;   // let GC reclaim
        if (twoScores == null) {
            Logging.info("  GPU kernel returned null after %d ms (infeasible)", gpuMs);
            return false;
        }
        Logging.info("  GPU kernel returned in %d ms", gpuMs);
        unpackTwoScores(twoScores, scoreArray, scoreArrayD, scoreArrayI, numSplits);
        return true;
    }

    /**
     * Convert the raw per-split doubled-score transport array from the GPU into
     * final per-split scores (score = 2·score / 2).
     *
     * <ul>
     *   <li>LONG  — each {@code twoScores[i]} is the exact integer 2·score.</li>
     *   <li>DOUBLE— the kernel stored the IEEE-754 bit pattern of the 2·score in
     *       the long slot ({@code __double_as_longlong}); recover with
     *       {@link Double#longBitsToDouble}.</li>
     *   <li>INT128— two longs per split: {@code [2i]} = low (unsigned),
     *       {@code [2i+1]} = high (signed).</li>
     * </ul>
     */
    private void unpackTwoScores(long[] twoScores, long[] scoreArray,
                                 double[] scoreArrayD, Int128[] scoreArrayI, int numSplits) {
        if (useInt128) {
            for (int i = 0; i < numSplits; i++) {
                long lo = twoScores[2 * i];
                long hi = twoScores[2 * i + 1];
                scoreArrayI[i] = new Int128(hi, lo).halve();   // 2·score → score
            }
        } else if (useDouble) {
            for (int i = 0; i < numSplits; i++)
                scoreArrayD[i] = Double.longBitsToDouble(twoScores[i]) / 2.0;
        } else {
            for (int i = 0; i < numSplits; i++)
                scoreArray[i] = twoScores[i] / 2L;
        }
    }

    /**
     * Legacy GPU path: smaller-side traversal, no prefix sums.
     *
     * Packs deduplicated binary rooted child partitions into the 8-int "parts" layout and calls
     * the one-thread-per-split kernel that counts each intersection by walking the
     * smaller range.  Uses the same splits and orderings/invIndex layout as the
     * prefix-sum path; builds NO node CSR and NO prefix working memory.
     */
    private boolean computeScoresGPUSmallerSide(List<BipartitionSplit> splitList,
                                                 PartitionTable partTable,
                                                 ClusterTable clusterTable,
                                                 List<Tree> clusterTrees,
                                                 List<Tree> partTrees,
                                                 int numGpuTrees,
                                                 long[] scoreArray,
                                                 double[] scoreArrayD,
                                                 Int128[] scoreArrayI,
                                                 int batchSizeHint,
                                                 double vramFraction) {
        int numSplits      = splitList.size();
        boolean splitTrees = (clusterTrees != partTrees);
        int partTreeOffset = splitTrees ? clusterTrees.size() : 0;

        int[] splitsData = buildSplitsData(splitList, clusterTable);
        int[][] rng      = buildSplitRangeData(splitList, clusterTable);
        int[] splitRangeMeta = rng[0], rangeData = rng[1];

        // --- parts (binary d==3): numParts * 8 ints; poly (d>3): separate CSR ---
        // treeIdx stored as (p.treeIndex + partTreeOffset) so the kernel reads the
        // original-tree half of the combined orderings/invIndex.
        // [treeIdx, lo1, hi1, lo2, hi2, sz1, sz2, frequency]
        List<PartitionTable.Entry> binEntries = new ArrayList<>();
        List<PartitionTable.Entry> polyEntries = new ArrayList<>();
        for (PartitionTable.Entry pe : partTable.entries()) {
            if (pe.exemplar.d == 3) binEntries.add(pe);
            else                    polyEntries.add(pe);
        }
        int numParts = binEntries.size(), numPolyParts = polyEntries.size();
        long polyBoundsLen = 0;
        for (PartitionTable.Entry pe : polyEntries) polyBoundsLen += pe.exemplar.d;
        int[] partsData = new int[numParts * 8];
        int[] ssPolyMeta        = new int[numPolyParts * 2];
        int[] ssPolyBoundOffset = new int[numPolyParts + 1];
        int[] ssPolyBounds      = new int[(int) polyBoundsLen];
        int j = 0;
        for (PartitionTable.Entry pe : binEntries) {
            Partition p = pe.exemplar;
            int base = j * 8;
            partsData[base + 0] = p.treeIndex + partTreeOffset;
            partsData[base + 1] = p.leftStart;
            partsData[base + 2] = p.leftEnd;
            partsData[base + 3] = p.rightStart;
            partsData[base + 4] = p.rightEnd;
            partsData[base + 5] = p.size1;
            partsData[base + 6] = p.size2;
            partsData[base + 7] = pe.frequency;
            j++;
        }
        int pj = 0, boundCur = 0;
        for (PartitionTable.Entry pe : polyEntries) {
            Partition p = pe.exemplar;
            ssPolyMeta[pj * 2]     = p.treeIndex + partTreeOffset;
            ssPolyMeta[pj * 2 + 1] = pe.frequency;
            ssPolyBoundOffset[pj] = boundCur;
            int k = p.d - 1;
            for (int i = 0; i < k; i++) ssPolyBounds[boundCur + i] = p.partStarts[i];
            ssPolyBounds[boundCur + k] = p.partEnds[k - 1];
            boundCur += p.d;
            pj++;
        }
        ssPolyBoundOffset[numPolyParts] = boundCur;

        int[][] oi      = buildOrderingsInvIndex(clusterTrees, partTrees, numGpuTrees);
        int[] orderings = oi[0], invIndex = oi[1];

        long t1 = System.nanoTime();
        long[] twoScores = GPUWeightCalculator.computeWeightsSmallerSideGPU(
            splitsData, splitRangeMeta, rangeData, partsData,
            ssPolyMeta, ssPolyBoundOffset, ssPolyBounds, orderings, invIndex,
            numSplits, numParts, numPolyParts, numGpuTrees, n, n,
            batchSizeHint, vramFraction, nativeScoreMode(),
            Config.getInstance().getGpuProgressIntervalSec());
        long gpuMs = (System.nanoTime() - t1) / 1_000_000;

        splitsData = null; partsData = null; orderings = null; invIndex = null;   // let GC reclaim
        if (twoScores == null) {
            Logging.info("  GPU kernel returned null after %d ms (infeasible)", gpuMs);
            return false;
        }
        Logging.info("  GPU kernel returned in %d ms", gpuMs);
        unpackTwoScores(twoScores, scoreArray, scoreArrayD, scoreArrayI, numSplits);
        return true;
    }

    // --- shared GPU input packers (identical layout for both kernels) ---

    /**
     * splits: numSplits * 10 ints.  Cluster treeIndex values are 0..k-1 (completed
     * trees, used for membership).
     * [aTree, aLo, aHi, aComp, aSize, bTree, bLo, bHi, bComp, bSize]
     */
    private int[] buildSplitsData(List<BipartitionSplit> splitList, ClusterTable clusterTable) {
        int numSplits = splitList.size();
        int[] splitsData = new int[numSplits * 10];
        for (int i = 0; i < numSplits; i++) {
            BipartitionSplit split = splitList.get(i);
            ClusterTable.Entry eA = clusterTable.get(split.lo);
            ClusterTable.Entry eB = clusterTable.get(split.hi);
            int base = i * 10;
            // Pack both single- and multi-range clusters. For a multi-range side the
            // kernel ignores [lo,hi) (it reads the split's range descriptor instead),
            // but tree/comp/size are still needed — left/right hold the bounding span.
            if (eA != null && eB != null) {
                Cluster cA = eA.exemplar, cB = eB.exemplar;
                splitsData[base + 0] = cA.treeIndex;
                splitsData[base + 1] = cA.left;
                splitsData[base + 2] = cA.right;
                splitsData[base + 3] = cA.complement ? 1 : 0;
                splitsData[base + 4] = cA.size;
                splitsData[base + 5] = cB.treeIndex;
                splitsData[base + 6] = cB.left;
                splitsData[base + 7] = cB.right;
                splitsData[base + 8] = cB.complement ? 1 : 0;
                splitsData[base + 9] = cB.size;
            }
            // else: all zeros → empty clusters → kernel yields score 0 for this split.
        }
        return splitsData;
    }

    /**
     * Per-split range descriptor + resident flat range array for the GPU two-tier
     * multi-range path.
     *   meta[i*4 + {0,1,2,3}] = {aRngOff, aRngCnt, bRngOff, bRngCnt}  (offsets in PAIRS)
     *   rangeData             = concatenated [lo,hi] pairs of every multi-range split side
     * A single-range side has count 0 — the kernel then uses the split's [lo,hi) fast path.
     * For runs with no multi-range clusters, meta is all-zero and rangeData is empty
     * (so the kernel's per-leaf membership is byte-identical to before).
     *
     * @return int[2][] = {meta (numSplits*4), rangeData}
     */
    private int[][] buildSplitRangeData(List<BipartitionSplit> splitList, ClusterTable clusterTable) {
        int numSplits = splitList.size();
        int[] meta = new int[numSplits * 4];
        java.util.ArrayList<Integer> ranges = new java.util.ArrayList<>(); // flat lo,hi pairs
        for (int i = 0; i < numSplits; i++) {
            BipartitionSplit sp = splitList.get(i);
            ClusterTable.Entry eA = clusterTable.get(sp.lo);
            ClusterTable.Entry eB = clusterTable.get(sp.hi);
            if (eA != null && eA.exemplar.isMultiRange()) {
                Cluster c = eA.exemplar;
                meta[i * 4 + 0] = ranges.size() / 2;   // offset in pairs
                meta[i * 4 + 1] = c.los.length;
                for (int j = 0; j < c.los.length; j++) { ranges.add(c.los[j]); ranges.add(c.his[j]); }
            }
            if (eB != null && eB.exemplar.isMultiRange()) {
                Cluster c = eB.exemplar;
                meta[i * 4 + 2] = ranges.size() / 2;
                meta[i * 4 + 3] = c.los.length;
                for (int j = 0; j < c.los.length; j++) { ranges.add(c.los[j]); ranges.add(c.his[j]); }
            }
        }
        int[] rangeData = new int[ranges.size()];
        for (int k = 0; k < ranges.size(); k++) rangeData[k] = ranges.get(k);
        return new int[][]{ meta, rangeData };
    }

    /**
     * orderings + invIndex: numGpuTrees * n ints each.
     *   orderings[t*n + pos]   = postorderArray[pos]
     *   invIndex [t*n + taxon] = positionMap[taxon]  (-1 if absent)
     *
     * Layout when splitTrees (clusterTrees != partTrees):
     *   slots 0..C-1     from clusterTrees  — cluster membership   (C = clusterTrees.size())
     *   slots C..C+P-1   from partTrees      — gene-tree leaves     (P = partTrees.size())
     * Otherwise the single list fills slots 0..k-1.
     *
     * @return int[2][] = {orderings, invIndex}
     */
    private int[][] buildOrderingsInvIndex(List<Tree> clusterTrees, List<Tree> partTrees,
                                            int numGpuTrees) {
        int numClusterTrees = clusterTrees.size();
        boolean splitTrees  = (clusterTrees != partTrees);
        int[] orderings = new int[numGpuTrees * n];
        int[] invIndex  = new int[numGpuTrees * n];
        Arrays.fill(invIndex, -1);
        for (int t = 0; t < numClusterTrees; t++) {
            Tree tree = clusterTrees.get(t);
            int base = t * n;
            for (int pos = 0; pos < tree.leafCount; pos++) orderings[base + pos] = tree.postorderArray[pos];
            for (int taxon = 0; taxon < n; taxon++)        invIndex[base + taxon] = tree.positionMap[taxon];
        }
        if (splitTrees) {
            int numPartTrees = partTrees.size();
            for (int t = 0; t < numPartTrees; t++) {
                Tree tree = partTrees.get(t);
                int base = (numClusterTrees + t) * n;
                for (int pos = 0; pos < tree.leafCount; pos++) orderings[base + pos] = tree.postorderArray[pos];
                for (int taxon = 0; taxon < n; taxon++)        invIndex[base + taxon] = tree.positionMap[taxon];
            }
        }
        return new int[][]{ orderings, invIndex };
    }

    // -------------------------------------------------------------------------
    // Per-tree internal-node CSR (rooted child partitions as leaf intervals)
    // -------------------------------------------------------------------------

    /**
     * Compact, per-exemplar-tree representation of the DEDUPLICATED gene-tree
     * rooted child partitions. Each unique partition is stored once, as a
     * contiguous leaf interval (lo, mid, hi) of its exemplar tree (M1 = [lo,mid)
     * and M2 = [mid,hi)), together with its frequency (how many gene-tree nodes
     * realize it). Taxa elsewhere in the exemplar tree are not part of the
     * partition. Entries are bucketed by exemplar tree so the
     * GPU kernel can build each tree's leaf prefix sums once and score only that
     * tree's unique rooted child partitions.
     *
     * This recovers cross-tree dedup at O(L) working memory (one tree's prefix
     * live at a time).
     * Scoring is bit-identical to summing triplet weights over every node,
     * because Σ_nodes w ≡ Σ_unique frequency·w.
     */
    private static final class NodeCSR {
        int[] nodeData;       // numBinUnique * 3   [lo, mid, hi] of the exemplar (d==3 only)
        int[] nodeFreq;       // numBinUnique       frequency (occurrence count)
        int[] nodeOffset;     // numTrees + 1       CSR row pointers (binary nodes by exemplar tree)
        int[] partLeafCount;  // numTrees           leaf count L per tree
        int   maxLeafCount;   // max L over trees with ≥1 exemplar (shared-mem sizing)
        int   totalNodes;     // numBinUnique + numPolyUnique  (for logging)
        // Polytomy (d>3) CSR — empty when no polytomous partitions.
        int[] polyTreeOffset;   // numTrees + 1     poly nodes bucketed by exemplar tree
        int[] polyBoundOffset;  // numPoly + 1      range into polyBounds (length d) per poly node
        int[] polyBounds;       // Σ d              concatenated boundary lists b[0..d-1]
        int[] polyFreq;         // numPoly          occurrence count
    }

    /**
     * Build the deduplicated node CSR from the already-computed PartitionTable,
     * bucketing unique rooted child partitions by their exemplar tree index.
     */
    private static NodeCSR buildDedupNodeCSR(PartitionTable partTable,
                                              List<Tree> partTrees) {
        int numTrees = partTrees.size();

        // Pass 1: count BINARY (d==3) and POLY (d>3) unique partitions per exemplar
        // tree, and the total poly boundary length.
        int[] nodeOffset     = new int[numTrees + 1];   // binary nodes per tree
        int[] polyTreeOffset = new int[numTrees + 1];   // poly nodes per tree
        long  polyBoundsLen  = 0;
        for (PartitionTable.Entry e : partTable.entries()) {
            if (e.exemplar.d == 3) nodeOffset[e.exemplar.treeIndex + 1]++;
            else { polyTreeOffset[e.exemplar.treeIndex + 1]++; polyBoundsLen += e.exemplar.d; }
        }
        for (int g = 0; g < numTrees; g++) {
            nodeOffset[g + 1]     += nodeOffset[g];
            polyTreeOffset[g + 1] += polyTreeOffset[g];
        }
        int total     = nodeOffset[numTrees];           // # binary unique
        int numPoly   = polyTreeOffset[numTrees];       // # poly unique
        if ((long) total * 3 > Integer.MAX_VALUE || polyBoundsLen > Integer.MAX_VALUE) {
            throw new IllegalStateException("Too many partitions for a single int[] CSR");
        }

        // Pass 2: scatter into binary nodeData and poly bound CSR, bucketed by tree.
        int[] nodeData = new int[total * 3];
        int[] nodeFreq = new int[total];
        int[] polyFreq        = new int[numPoly];
        int[] polyBoundOffset = new int[numPoly + 1];
        int[] polyBounds      = new int[(int) polyBoundsLen];
        int[] binCursor       = nodeOffset.clone();     // per-tree binary write cursor
        int[] polyCursor      = polyTreeOffset.clone(); // per-tree poly write cursor

        // Pre-fill polyBoundOffset by walking poly nodes in the SAME scatter order.
        // We fill it incrementally during scatter via a running bound cursor per slot.
        int[] polyDeg = new int[numPoly];               // degree per poly slot (for offsets)
        for (PartitionTable.Entry e : partTable.entries()) {
            Partition p = e.exemplar;
            if (p.d == 3) {
                int pos = binCursor[p.treeIndex]++;
                int b   = pos * 3;
                nodeData[b]     = p.leftStart;          // lo
                nodeData[b + 1] = p.leftEnd;            // mid  (= p.rightStart)
                nodeData[b + 2] = p.rightEnd;           // hi
                nodeFreq[pos]   = e.frequency;
            } else {
                int pos = polyCursor[p.treeIndex]++;
                polyFreq[pos] = e.frequency;
                polyDeg[pos]  = p.d;
            }
        }
        // Build polyBoundOffset (prefix sum of degrees) then scatter boundary lists.
        for (int pn = 0; pn < numPoly; pn++) polyBoundOffset[pn + 1] = polyBoundOffset[pn] + polyDeg[pn];
        int[] polyCursor2 = polyTreeOffset.clone();
        for (PartitionTable.Entry e : partTable.entries()) {
            Partition p = e.exemplar;
            if (p.d == 3) continue;
            int pos  = polyCursor2[p.treeIndex]++;
            int base = polyBoundOffset[pos];
            int k    = p.d - 1;                          // # child intervals
            for (int i = 0; i < k; i++) polyBounds[base + i] = p.partStarts[i];
            polyBounds[base + k] = p.partEnds[k - 1];    // final boundary = hi
        }

        // Per-tree leaf counts; maxLeaf over trees that need a prefix (binary OR poly).
        int[] partLeafCount = new int[numTrees];
        int   maxLeaf = 0;
        for (int g = 0; g < numTrees; g++) {
            int L = partTrees.get(g).leafCount;
            partLeafCount[g] = L;
            boolean hasWork = nodeOffset[g + 1] > nodeOffset[g] || polyTreeOffset[g + 1] > polyTreeOffset[g];
            if (hasWork && L > maxLeaf) maxLeaf = L;
        }

        NodeCSR csr = new NodeCSR();
        csr.nodeData        = nodeData;
        csr.nodeFreq        = nodeFreq;
        csr.nodeOffset      = nodeOffset;
        csr.partLeafCount   = partLeafCount;
        csr.maxLeafCount    = maxLeaf;
        csr.totalNodes      = total + numPoly;
        csr.polyTreeOffset  = polyTreeOffset;
        csr.polyBoundOffset = polyBoundOffset;
        csr.polyBounds      = polyBounds;
        csr.polyFreq        = polyFreq;
        return csr;
    }

    // -------------------------------------------------------------------------
    // BITSET path (low-taxa fast option; CPU + GPU).
    //
    // Every cluster and every gene-tree part is materialized ONCE as a global-taxon
    // bitset of W = ceil(n/64) 64-bit words.  Each core intersection |X ∩ Y| is then
    // popcount(X & Y) over W words — O(1) for small n, and independent of which tree
    // either set came from (both are keyed by global taxon id, so there is no
    // cross-tree coordinate walk and no orderings/invIndex in the score loop).  The
    // intersection derivation and rooted-triplet math use the exact same helpers the other
    // methods use, so scores are bit-identical.  Best when n is small and the gene
    // count is large; the per-intersection cost grows with W as n grows.
    // -------------------------------------------------------------------------

    /**
     * Precomputed global-taxon bitsets for the whole weight computation.
     *
     * A "cluster pool" holds one deduplicated bitset per distinct candidate-cluster
     * (cid 0 is a reserved all-zero, size-0 empty cluster used for splits whose side
     * is absent from the ClusterTable — those score 0, matching the range paths).
     * Parts are split into binary (M1,M2 bitsets) and polytomous (a CSR of child
     * bitsets).  No outside-node or present-taxa bitset is needed: rooted triplet
     * support is a function only of the actual children of the gene-tree node.
     */
    private static final class BitsetData {
        int      W;                 // 64-bit words per set = ceil(n/64)
        int      numClusters;
        long[]   clusterBits;       // numClusters * W
        int[]    splitCid;          // numSplits * 4  [aCid, bCid, aSize, bSize]
        int      numBin;
        long[]   partM1, partM2;    // numBin * W  each
        int[]    partFreq;          // numBin
        int      numPoly;
        int[]    polyMeta;          // numPoly * 2  [childCount, freq]
        int[]    polyChildOffset;   // numPoly + 1  CSR row pointers
        long[]   polyChildBits;     // totalChildren * W
    }

    /** Set bits for the taxa in postorder range [lo,hi) of tree t into arr[base..base+W). */
    private static void setRangeBits(long[] arr, int base, Tree t, int lo, int hi) {
        int[] po = t.postorderArray;
        for (int pos = lo; pos < hi; pos++) {
            int tax = po[pos];
            arr[base + (tax >>> 6)] |= 1L << (tax & 63);
        }
    }

    /** Write one cluster's global-taxon bitset into arr[base..base+W) (complement-aware). */
    private void buildClusterBitsInto(long[] arr, int base, Cluster c, List<Tree> clusterTrees, int W) {
        Tree t = clusterTrees.get(c.treeIndex);
        if (c.isMultiRange()) {
            for (int j = 0; j < c.los.length; j++) setRangeBits(arr, base, t, c.los[j], c.his[j]);
        } else {
            setRangeBits(arr, base, t, c.left, c.right);
        }
        if (c.complement) {                       // A = S \ range, complement over ALL n taxa
            for (int k = 0; k < W; k++) arr[base + k] = ~arr[base + k];
            int rem = n & 63;
            if (rem != 0) arr[base + W - 1] &= (1L << rem) - 1;   // clear bits ≥ n in the top word
        }
    }

    /**
     * Materialize all cluster/part/Lg bitsets once; shared by the CPU and GPU bitset paths.
     *
     * Each bitset occupies a disjoint slice of its output array, so the three fills are
     * done in parallel across threads.  Only the cluster-id assignment (a HashMap dedup)
     * is sequential; the actual cluster bitset construction is then parallel.
     */
    private BitsetData buildBitsetData(List<BipartitionSplit> splitList,
                                       PartitionTable partTable, ClusterTable clusterTable,
                                       List<Tree> clusterTrees, List<Tree> partTrees) {
        int W = (n + 63) >>> 6;
        int numSplits = splitList.size();

        // --- Cluster pool: phase 1 sequential id assignment (dedup by ClusterHash) ---
        // cid 0 = empty cluster (all-zero, size 0): splits whose side is absent score 0.
        Map<ClusterHash, Integer> cidMap = new HashMap<>();
        List<Cluster> clusterExemplars = new ArrayList<>();
        clusterExemplars.add(null);               // cid 0 → empty
        int[] splitCid = new int[numSplits * 4];
        ClusterHash[] side = new ClusterHash[2];
        for (int i = 0; i < numSplits; i++) {
            BipartitionSplit sp = splitList.get(i);
            side[0] = sp.lo; side[1] = sp.hi;
            for (int s = 0; s < 2; s++) {
                ClusterTable.Entry e = clusterTable.get(side[s]);
                int cid = 0, sz = 0;
                if (e != null) {
                    Integer id = cidMap.get(side[s]);
                    if (id == null) {
                        id = clusterExemplars.size();
                        clusterExemplars.add(e.exemplar);
                        cidMap.put(side[s], id);
                    }
                    cid = id; sz = e.exemplar.size;
                }
                splitCid[i * 4 + s]     = cid;
                splitCid[i * 4 + 2 + s] = sz;
            }
        }
        int numClusters = clusterExemplars.size();
        if ((long) numClusters * W > Integer.MAX_VALUE)
            throw new IllegalStateException("Bitset cluster pool too large for a single long[]");
        // Phase 2: parallel fill of the deduplicated cluster bitsets (disjoint slices).
        long[] clusterBits = new long[numClusters * W];
        Threading.processRangeParallel(numClusters, cid -> {
            if (cid == 0) return;                 // cid 0 stays all-zero
            buildClusterBitsInto(clusterBits, cid * W, clusterExemplars.get(cid), clusterTrees, W);
        });

        // --- Part bitsets: split entries into binary / poly, prefix the poly child cursor,
        //     then fill each group in parallel (each part owns a disjoint slice). ---
        List<PartitionTable.Entry> binEntries  = new ArrayList<>();
        List<PartitionTable.Entry> polyEntries = new ArrayList<>();
        for (PartitionTable.Entry e : partTable.entries()) {
            if (e.exemplar.d == 3) binEntries.add(e);
            else                   polyEntries.add(e);
        }
        int numBin = binEntries.size(), numPoly = polyEntries.size();
        int[] polyChildOffset = new int[numPoly + 1];
        for (int pi = 0; pi < numPoly; pi++)
            polyChildOffset[pi + 1] = polyChildOffset[pi] + (polyEntries.get(pi).exemplar.d - 1);
        long childTotal = polyChildOffset[numPoly];
        if ((long) numBin * W > Integer.MAX_VALUE || childTotal * W > Integer.MAX_VALUE)
            throw new IllegalStateException("Bitset part pool too large for a single long[]");

        long[] partM1 = new long[numBin * W];
        long[] partM2 = new long[numBin * W];
        int[]  partFreq = new int[numBin];
        int[]  polyMeta = new int[numPoly * 2];
        long[] polyChildBits = new long[(int) childTotal * W];

        Threading.processRangeParallel(numBin, bi -> {
            PartitionTable.Entry e = binEntries.get(bi);
            Partition p = e.exemplar;
            Tree tGT = partTrees.get(p.treeIndex);
            setRangeBits(partM1, bi * W, tGT, p.leftStart,  p.leftEnd);
            setRangeBits(partM2, bi * W, tGT, p.rightStart, p.rightEnd);
            partFreq[bi] = e.frequency;
        });
        Threading.processRangeParallel(numPoly, pi -> {
            PartitionTable.Entry e = polyEntries.get(pi);
            Partition p = e.exemplar;
            Tree tGT = partTrees.get(p.treeIndex);
            int d = p.d, childCur = polyChildOffset[pi];
            for (int c = 0; c < d - 1; c++) {
                setRangeBits(polyChildBits, (childCur + c) * W, tGT, p.partStarts[c], p.partEnds[c]);
            }
            int m = pi * 2;
            polyMeta[m]     = d - 1;
            polyMeta[m + 1] = e.frequency;
        });

        BitsetData bd = new BitsetData();
        bd.W = W;
        bd.numClusters = numClusters;
        bd.clusterBits = clusterBits;
        bd.splitCid = splitCid;
        bd.numBin = numBin;
        bd.partM1 = partM1;
        bd.partM2 = partM2;
        bd.partFreq = partFreq;
        bd.numPoly = numPoly;
        bd.polyMeta = polyMeta;
        bd.polyChildOffset = polyChildOffset;
        bd.polyChildBits = polyChildBits;
        return bd;
    }

    /** popcount( X[xBase..] & Y[yBase..] ) over W words. */
    private static int popAnd(long[] x, int xBase, long[] y, int yBase, int W) {
        int c = 0;
        for (int k = 0; k < W; k++) c += Long.bitCount(x[xBase + k] & y[yBase + k]);
        return c;
    }

    // --- CPU bitset scorers (mirror computeScore / D / I; reuse score helpers) ---

    private long computeScoreBitset(int i, BitsetData bd, int totalN) {
        int W = bd.W;
        int aCid = bd.splitCid[i * 4], bCid = bd.splitCid[i * 4 + 1];
        int aSize = bd.splitCid[i * 4 + 2], bSize = bd.splitCid[i * 4 + 3];
        if (totalN - aSize - bSize < 0) return 0L;
        int aB = aCid * W, bB = bCid * W;
        long twoScore = 0L;

        for (int j = 0; j < bd.numBin; j++) {
            int mW = j * W;
            int a0 = popAnd(bd.clusterBits, aB, bd.partM1, mW, W);
            int a1 = popAnd(bd.clusterBits, aB, bd.partM2, mW, W);
            int b0 = popAnd(bd.clusterBits, bB, bd.partM1, mW, W);
            int b1 = popAnd(bd.clusterBits, bB, bd.partM2, mW, W);
            twoScore += (long) bd.partFreq[j] * computeTwoQI(a0, a1, b0, b1);
        }
        for (int pn = 0; pn < bd.numPoly; pn++) {
            int m = pn * 2, childCount = bd.polyMeta[m], freq = bd.polyMeta[m + 1];
            int[][] parts = polyChildIntersectionsBitset(bd, pn, childCount, aB, bB);
            twoScore += (long) freq * polyTwoQILong(parts[0], parts[1], childCount);
        }
        return twoScore / 2;
    }

    private double computeScoreBitsetD(int i, BitsetData bd, int totalN) {
        int W = bd.W;
        int aCid = bd.splitCid[i * 4], bCid = bd.splitCid[i * 4 + 1];
        int aSize = bd.splitCid[i * 4 + 2], bSize = bd.splitCid[i * 4 + 3];
        if (totalN - aSize - bSize < 0) return 0.0;
        int aB = aCid * W, bB = bCid * W;
        double twoScore = 0.0;

        for (int j = 0; j < bd.numBin; j++) {
            int mW = j * W;
            int a0 = popAnd(bd.clusterBits, aB, bd.partM1, mW, W);
            int a1 = popAnd(bd.clusterBits, aB, bd.partM2, mW, W);
            int b0 = popAnd(bd.clusterBits, bB, bd.partM1, mW, W);
            int b1 = popAnd(bd.clusterBits, bB, bd.partM2, mW, W);
            twoScore += (double) bd.partFreq[j] * computeTwoQIDouble(a0, a1, b0, b1);
        }
        for (int pn = 0; pn < bd.numPoly; pn++) {
            int m = pn * 2, childCount = bd.polyMeta[m], freq = bd.polyMeta[m + 1];
            int[][] parts = polyChildIntersectionsBitset(bd, pn, childCount, aB, bB);
            twoScore += (double) freq * polyTwoQIDouble(parts[0], parts[1], childCount);
        }
        return twoScore / 2.0;
    }

    private Int128 computeScoreBitsetI(int i, BitsetData bd, int totalN) {
        int W = bd.W;
        int aCid = bd.splitCid[i * 4], bCid = bd.splitCid[i * 4 + 1];
        int aSize = bd.splitCid[i * 4 + 2], bSize = bd.splitCid[i * 4 + 3];
        if (totalN - aSize - bSize < 0) return Int128.ZERO;
        int aB = aCid * W, bB = bCid * W;
        Int128 twoScore = Int128.ZERO;

        for (int j = 0; j < bd.numBin; j++) {
            int mW = j * W;
            int a0 = popAnd(bd.clusterBits, aB, bd.partM1, mW, W);
            int a1 = popAnd(bd.clusterBits, aB, bd.partM2, mW, W);
            int b0 = popAnd(bd.clusterBits, bB, bd.partM1, mW, W);
            int b1 = popAnd(bd.clusterBits, bB, bd.partM2, mW, W);
            twoScore = twoScore.add(computeTwoQIInt128(a0, a1, b0, b1).mulScalar(bd.partFreq[j]));
        }
        for (int pn = 0; pn < bd.numPoly; pn++) {
            int m = pn * 2, childCount = bd.polyMeta[m], freq = bd.polyMeta[m + 1];
            int[][] parts = polyChildIntersectionsBitset(bd, pn, childCount, aB, bB);
            twoScore = twoScore.add(polyTwoQIInt128(parts[0], parts[1], childCount).mulScalar(freq));
        }
        return twoScore.halve();
    }

    /**
     * Bitset analogue of {@link #polyChildIntersections}: build the two child-only
     * intersection rows for a polytomous part.
     */
    private int[][] polyChildIntersectionsBitset(BitsetData bd, int pn, int childCount,
                                                  int aB, int bB) {
        int W = bd.W;
        int cbeg = bd.polyChildOffset[pn];
        int[] a = new int[childCount], b = new int[childCount];
        for (int ci = 0; ci < childCount; ci++) {
            int mW = (cbeg + ci) * W;
            a[ci] = popAnd(bd.clusterBits, aB, bd.polyChildBits, mW, W);
            b[ci] = popAnd(bd.clusterBits, bB, bd.polyChildBits, mW, W);
        }
        return new int[][]{ a, b };
    }

    /** CPU bitset path: build bitsets once, then score splits in parallel. */
    private void computeScoresCPUBitset(List<BipartitionSplit> splitList,
                                        PartitionTable partTable, ClusterTable clusterTable,
                                        List<Tree> clusterTrees, List<Tree> partTrees,
                                        long[] scoreArray, double[] scoreArrayD, Int128[] scoreArrayI) {
        BitsetData bd = buildBitsetData(splitList, partTable, clusterTable, clusterTrees, partTrees);
        int numSplits = splitList.size();
        int totalN = n;
        Logging.info("Weight table: CPU path (bitset)  splits=%d  W=%d  clusters=%d  bin=%d  poly=%d",
            numSplits, bd.W, bd.numClusters, bd.numBin, bd.numPoly);
        java.util.concurrent.atomic.AtomicInteger wDone = new java.util.concurrent.atomic.AtomicInteger(0);
        ProgressBar wBar = new ProgressBar("Scoring splits (CPU bitset)", numSplits);
        Threading.processRangeParallel(numSplits, idx -> {
            if (useInt128)      scoreArrayI[idx] = computeScoreBitsetI(idx, bd, totalN);
            else if (useDouble) scoreArrayD[idx] = computeScoreBitsetD(idx, bd, totalN);
            else                scoreArray[idx]  = computeScoreBitset(idx, bd, totalN);
            wBar.update(wDone.incrementAndGet());
        });
        wBar.done();
    }

    /** GPU bitset path: build bitsets once, resolve batch size, call the native kernel. */
    private boolean computeScoresGPUBitset(List<BipartitionSplit> splitList,
                                           PartitionTable partTable, ClusterTable clusterTable,
                                           List<Tree> clusterTrees, List<Tree> partTrees,
                                           long[] scoreArray, double[] scoreArrayD, Int128[] scoreArrayI) {
        Config cfg = Config.getInstance();
        BitsetData bd = buildBitsetData(splitList, partTable, clusterTable, clusterTrees, partTrees);
        int numSplits = splitList.size();

        // Resident device memory (for vram-control-factor sizing + logging only).
        long residentMem = 8L * (bd.clusterBits.length + bd.partM1.length + bd.partM2.length
                                 + bd.polyChildBits.length)
                         + 4L * (bd.partFreq.length + bd.polyMeta.length
                                 + bd.polyChildOffset.length);
        long perSplit = 4L * Integer.BYTES + (useInt128 ? 2L : 1L) * Long.BYTES;   // 4 ints in + score out

        // Batch-size hint — same priority as the other GPU paths.
        int batchSizeHint; String batchDesc;
        if (!cfg.isGpuBatch()) {
            batchSizeHint = -1; batchDesc = "off (single launch)";
        } else if (cfg.getGpuNumBatches() > 0) {
            int N = cfg.getGpuNumBatches();
            batchSizeHint = (numSplits + N - 1) / N;
            batchDesc = N + " batches → batchSize=" + batchSizeHint;
        } else if (cfg.getGpuBatchSize() > 0) {
            batchSizeHint = cfg.getGpuBatchSize();
            batchDesc = "explicit batchSize=" + batchSizeHint;
        } else if (cfg.isGpuVramControlFactorSet()) {
            double F = cfg.getGpuVramControlFactor();
            long batchMem = (long) (F * residentMem);
            batchSizeHint = (int) Math.max(1, Math.min(numSplits, batchMem / perSplit));
            batchDesc = String.format("vram-control-factor=%.3f  resident=%.1f MB  batch→%d",
                F, residentMem / 1e6, (numSplits + batchSizeHint - 1) / Math.max(1, batchSizeHint));
        } else {
            batchSizeHint = 0;
            batchDesc = String.format("auto (free-VRAM adaptive, occupancy=%.0f%%)",
                cfg.getGpuVramFraction() * 100);
        }

        Logging.info("Weight table: GPU path (bitset)  splits=%d  W=%d  clusters=%d  bin=%d  poly=%d  resident=%.1f MB  batching=%s",
            numSplits, bd.W, bd.numClusters, bd.numBin, bd.numPoly,
            residentMem / 1e6, batchDesc);

        long t1 = System.nanoTime();
        long[] twoScores = GPUWeightCalculator.computeWeightsBitsetGPU(
            bd.splitCid, bd.clusterBits, bd.partM1, bd.partM2, bd.partFreq,
            bd.polyMeta, bd.polyChildOffset, bd.polyChildBits,
            numSplits, bd.numClusters, bd.numBin, bd.numPoly, bd.W, n,
            batchSizeHint, cfg.getGpuVramFraction(), nativeScoreMode(),
            cfg.getGpuProgressIntervalSec());
        long gpuMs = (System.nanoTime() - t1) / 1_000_000;

        if (twoScores == null) {
            Logging.info("  GPU bitset kernel returned null after %d ms (infeasible)", gpuMs);
            return false;
        }
        Logging.info("  GPU bitset kernel returned in %d ms", gpuMs);
        unpackTwoScores(twoScores, scoreArray, scoreArrayD, scoreArrayI, numSplits);
        return true;
    }

    // -------------------------------------------------------------------------
    // SIMPLE-TREE-WALK path (many-candidate fast option; CPU + GPU).
    //
    // One thread per split walks a resident flat postorder token stream of all gene
    // trees sequentially, maintaining a small O(n) per-thread stack of
    // (nA,nB) = (|node∩A|, |node∩B|) pairs. Every rooted internal
    // node, including the gene-tree root, is scored from its actual children — the
    // same node set and the same rooted-triplet helpers as the other methods, with
    // no dedup, so raw Σ over nodes is bit-identical to deduped Σ frequency·weight.
    // -------------------------------------------------------------------------

    /** GPU per-thread postorder stack cap — mirrors WB_TW_STACK_CAP in stelarx_weight.cu. */
    private static final int TW_GPU_STACK_CAP = 512;

    private static final class TreeWalkData {
        int      W;
        int      numClusters;
        long[]   clusterBits;     // numClusters * W
        int[]    splitCid;        // numSplits * 4  [aCid, bCid, aSize, bSize]
        int      numTrees;
        int[]    nodeStream;      // flat postorder tokens (leaf=taxon≥0, internal=-childCount)
        int[]    treeNodeOffset;  // numTrees + 1  CSR row pointers
        int      maxFrontier;     // exact max postorder stack entries over all trees
        int      frontierTree;    // tree index attaining maxFrontier
        int      frontierLeaves;  // leaf count of frontierTree
    }

    /** Exact tree-walk stack requirement for the current child ordering. */
    private static final class TreeWalkFrontier {
        int entries = 1;
        int treeIndex = -1;
        int leafCount = 0;
    }

    /**
     * Maximum evaluation-stack size of this subtree's postorder token sequence.
     * Each completed child leaves one value on the stack; the internal-node token
     * then reduces all child values to one without increasing the peak.
     */
    private static int treeWalkFrontier(TreeNode node) {
        if (node.isLeaf()) return 1;
        int held = 0;
        int peak = 0;
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) {
                peak = Math.max(peak, held + treeWalkFrontier(child));
                held++;
            }
        } else {
            peak = Math.max(peak, treeWalkFrontier(node.left));
            held = 1;
            peak = Math.max(peak, held + treeWalkFrontier(node.right));
        }
        return Math.max(1, peak);
    }

    /** Measure the exact maximum frontier across all scoring gene trees. */
    private static TreeWalkFrontier measureTreeWalkFrontier(List<Tree> trees) {
        TreeWalkFrontier result = new TreeWalkFrontier();
        for (int g = 0; g < trees.size(); g++) {
            Tree tree = trees.get(g);
            int frontier = treeWalkFrontier(tree.root);
            if (frontier > result.entries || result.treeIndex < 0) {
                result.entries = frontier;
                result.treeIndex = g;
                result.leafCount = tree.leafCount;
            }
        }
        return result;
    }

    /** Independently replay emitted tokens to audit the structural measurement. */
    private static int tokenStreamFrontier(int[] stream, int begin, int end) {
        int top = 0;
        int peak = 0;
        for (int i = begin; i < end; i++) {
            int token = stream[i];
            if (token >= 0) {
                top++;
                if (top > peak) peak = top;
            } else {
                int childCount = -token;
                if (childCount > top) {
                    throw new IllegalStateException("Malformed tree-walk token stream: childCount="
                        + childCount + " exceeds stack top=" + top);
                }
                top = top - childCount + 1;
            }
        }
        return Math.max(1, peak);
    }

    /** Count postorder tokens: every leaf and every rooted internal node, root included. */
    private int countTreeTokens(TreeNode node) {
        if (node.isLeaf()) return 1;
        int c = 0;
        if (node.isPolytomous()) for (TreeNode ch : node.children) c += countTreeTokens(ch);
        else { c += countTreeTokens(node.left); c += countTreeTokens(node.right); }
        c += 1;
        return c;
    }

    /** Emit postorder tokens; leaf → taxonId (≥0), internal → -childCount. */
    private int fillTreeTokens(TreeNode node, int[] stream, int cursor) {
        if (node.isLeaf()) { stream[cursor++] = node.taxonId; return cursor; }
        int k;
        if (node.isPolytomous()) {
            k = node.children.length;
            for (TreeNode ch : node.children) cursor = fillTreeTokens(ch, stream, cursor);
        } else {
            k = 2;
            cursor = fillTreeTokens(node.left,  stream, cursor);
            cursor = fillTreeTokens(node.right, stream, cursor);
        }
        stream[cursor++] = -k;
        return cursor;
    }

    /** Materialize the cluster pool and postorder token stream (parallel fills). */
    private TreeWalkData buildTreeWalkData(List<BipartitionSplit> splitList, ClusterTable clusterTable,
                                           List<Tree> clusterTrees, List<Tree> partTrees,
                                           TreeWalkFrontier frontier) {
        int W = (n + 63) >>> 6;
        int numSplits = splitList.size();
        int numTrees = partTrees.size();

        // --- Cluster pool (dedup by ClusterHash); cid 0 = empty (all-zero, size 0) ---
        Map<ClusterHash, Integer> cidMap = new HashMap<>();
        List<Cluster> clusterExemplars = new ArrayList<>();
        clusterExemplars.add(null);               // cid 0
        int[] splitCid = new int[numSplits * 4];
        ClusterHash[] side = new ClusterHash[2];
        for (int i = 0; i < numSplits; i++) {
            BipartitionSplit sp = splitList.get(i);
            side[0] = sp.lo; side[1] = sp.hi;
            for (int s = 0; s < 2; s++) {
                ClusterTable.Entry e = clusterTable.get(side[s]);
                int cid = 0, sz = 0;
                if (e != null) {
                    Integer id = cidMap.get(side[s]);
                    if (id == null) {
                        id = clusterExemplars.size();
                        clusterExemplars.add(e.exemplar);
                        cidMap.put(side[s], id);
                    }
                    cid = id; sz = e.exemplar.size;
                }
                splitCid[i * 4 + s]     = cid;
                splitCid[i * 4 + 2 + s] = sz;
            }
        }
        int numClusters = clusterExemplars.size();
        if ((long) numClusters * W > Integer.MAX_VALUE)
            throw new IllegalStateException("Tree-walk cluster pool too large for a single long[]");
        long[] clusterBits = new long[numClusters * W];
        Threading.processRangeParallel(numClusters, cid -> {
            if (cid == 0) return;
            buildClusterBitsInto(clusterBits, cid * W, clusterExemplars.get(cid), clusterTrees, W);
        });

        // --- Postorder token stream: count → offsets → parallel fill (disjoint segments) ---
        int[] treeNodeOffset = new int[numTrees + 1];
        for (int g = 0; g < numTrees; g++) {
            Tree t = partTrees.get(g);
            long tok = (long) treeNodeOffset[g] + countTreeTokens(t.root);
            if (tok > Integer.MAX_VALUE)
                throw new IllegalStateException("Tree-walk node stream too large for a single int[]");
            treeNodeOffset[g + 1] = (int) tok;
        }
        int[] nodeStream = new int[treeNodeOffset[numTrees]];
        Threading.processRangeParallel(numTrees, g ->
            fillTreeTokens(partTrees.get(g).root, nodeStream, treeNodeOffset[g]));

        int replayMax = 1;
        int replayTree = -1;
        for (int g = 0; g < numTrees; g++) {
            int measured = tokenStreamFrontier(nodeStream, treeNodeOffset[g], treeNodeOffset[g + 1]);
            if (measured > replayMax || replayTree < 0) {
                replayMax = measured;
                replayTree = g;
            }
        }
        if (replayMax != frontier.entries) {
            throw new IllegalStateException("Tree-walk frontier mismatch: structural="
                + frontier.entries + " tokenReplay=" + replayMax + " (tree=" + replayTree + ")");
        }

        TreeWalkData d = new TreeWalkData();
        d.W = W; d.numClusters = numClusters; d.clusterBits = clusterBits; d.splitCid = splitCid;
        d.numTrees = numTrees; d.nodeStream = nodeStream; d.treeNodeOffset = treeNodeOffset;
        d.maxFrontier = frontier.entries; d.frontierTree = frontier.treeIndex;
        d.frontierLeaves = frontier.leafCount;
        return d;
    }

    // --- CPU tree-walk scorers (mirror computeScore / D / I; reuse score helpers) ---

    private long computeScoreTreeWalk(int idx, TreeWalkData d, int totalN, int[] stack) {
        int W = d.W;
        int aCid = d.splitCid[idx * 4], bCid = d.splitCid[idx * 4 + 1];
        int aSize = d.splitCid[idx * 4 + 2], bSize = d.splitCid[idx * 4 + 3];
        if (totalN - aSize - bSize < 0) return 0L;
        int aB = aCid * W, bB = bCid * W;
        long[] cb = d.clusterBits;
        int[] ns = d.nodeStream, off = d.treeNodeOffset;
        long twoScore = 0L;

        for (int g = 0; g < d.numTrees; g++) {
            int segBeg = off[g], segEnd = off[g + 1];
            if (segBeg == segEnd) continue;

            int top = 0;
            for (int i = segBeg; i < segEnd; i++) {
                int tok = ns[i];
                if (tok >= 0) {
                    int inA = (int) ((cb[aB + (tok >>> 6)] >>> (tok & 63)) & 1L);
                    int inB = (int) ((cb[bB + (tok >>> 6)] >>> (tok & 63)) & 1L);
                    int e = top * 2; stack[e] = inA; stack[e + 1] = inB; top++;
                } else {
                    int k = -tok, cbase = top - k;
                    if (k == 2) {
                        int e0 = cbase * 2;
                        int a0 = stack[e0], b0 = stack[e0 + 1];
                        int a1 = stack[e0 + 2], b1 = stack[e0 + 3];
                        twoScore += computeTwoQI(a0, a1, b0, b1);
                        stack[e0] = a0 + a1; stack[e0 + 1] = b0 + b1;
                        top = cbase + 1;
                    } else {
                        long sumA = 0, sumB = 0;
                        for (int j = 0; j < k; j++) {
                            int e = (cbase + j) * 2;
                            sumA += stack[e]; sumB += stack[e + 1];
                        }
                        for (int j = 0; j < k; j++) {
                            int e = (cbase + j) * 2;
                            long aj = stack[e], bj = stack[e + 1];
                            twoScore += aj * (aj - 1L) * (sumB - bj)
                                      + bj * (bj - 1L) * (sumA - aj);
                        }
                        stack[cbase * 2] = (int) sumA; stack[cbase * 2 + 1] = (int) sumB;
                        top = cbase + 1;
                    }
                }
            }
        }
        return twoScore / 2;
    }

    private double computeScoreTreeWalkD(int idx, TreeWalkData d, int totalN, int[] stack) {
        int W = d.W;
        int aCid = d.splitCid[idx * 4], bCid = d.splitCid[idx * 4 + 1];
        int aSize = d.splitCid[idx * 4 + 2], bSize = d.splitCid[idx * 4 + 3];
        if (totalN - aSize - bSize < 0) return 0.0;
        int aB = aCid * W, bB = bCid * W;
        long[] cb = d.clusterBits;
        int[] ns = d.nodeStream, off = d.treeNodeOffset;
        double twoScore = 0.0;

        for (int g = 0; g < d.numTrees; g++) {
            int segBeg = off[g], segEnd = off[g + 1];
            if (segBeg == segEnd) continue;
            int top = 0;
            for (int i = segBeg; i < segEnd; i++) {
                int tok = ns[i];
                if (tok >= 0) {
                    int inA = (int) ((cb[aB + (tok >>> 6)] >>> (tok & 63)) & 1L);
                    int inB = (int) ((cb[bB + (tok >>> 6)] >>> (tok & 63)) & 1L);
                    int e = top * 2; stack[e] = inA; stack[e + 1] = inB; top++;
                } else {
                    int k = -tok, cbase = top - k;
                    if (k == 2) {
                        int e0 = cbase * 2;
                        int a0 = stack[e0], b0 = stack[e0 + 1];
                        int a1 = stack[e0 + 2], b1 = stack[e0 + 3];
                        twoScore += computeTwoQIDouble(a0, a1, b0, b1);
                        stack[e0] = a0 + a1; stack[e0 + 1] = b0 + b1;
                        top = cbase + 1;
                    } else {
                        long sumA = 0, sumB = 0;
                        for (int j = 0; j < k; j++) {
                            int e = (cbase + j) * 2;
                            sumA += stack[e]; sumB += stack[e + 1];
                        }
                        for (int j = 0; j < k; j++) {
                            int e = (cbase + j) * 2;
                            double aj = stack[e], bj = stack[e + 1];
                            twoScore += aj * (aj - 1.0) * (sumB - bj)
                                      + bj * (bj - 1.0) * (sumA - aj);
                        }
                        stack[cbase * 2] = (int) sumA; stack[cbase * 2 + 1] = (int) sumB;
                        top = cbase + 1;
                    }
                }
            }
        }
        return twoScore / 2.0;
    }

    private Int128 computeScoreTreeWalkI(int idx, TreeWalkData d, int totalN, int[] stack) {
        int W = d.W;
        int aCid = d.splitCid[idx * 4], bCid = d.splitCid[idx * 4 + 1];
        int aSize = d.splitCid[idx * 4 + 2], bSize = d.splitCid[idx * 4 + 3];
        if (totalN - aSize - bSize < 0) return Int128.ZERO;
        int aB = aCid * W, bB = bCid * W;
        long[] cb = d.clusterBits;
        int[] ns = d.nodeStream, off = d.treeNodeOffset;
        Int128 twoScore = Int128.ZERO;

        for (int g = 0; g < d.numTrees; g++) {
            int segBeg = off[g], segEnd = off[g + 1];
            if (segBeg == segEnd) continue;
            int top = 0;
            for (int i = segBeg; i < segEnd; i++) {
                int tok = ns[i];
                if (tok >= 0) {
                    int inA = (int) ((cb[aB + (tok >>> 6)] >>> (tok & 63)) & 1L);
                    int inB = (int) ((cb[bB + (tok >>> 6)] >>> (tok & 63)) & 1L);
                    int e = top * 2; stack[e] = inA; stack[e + 1] = inB; top++;
                } else {
                    int k = -tok, cbase = top - k;
                    if (k == 2) {
                        int e0 = cbase * 2;
                        int a0 = stack[e0], b0 = stack[e0 + 1];
                        int a1 = stack[e0 + 2], b1 = stack[e0 + 3];
                        twoScore = twoScore.add(computeTwoQIInt128(a0, a1, b0, b1));
                        stack[e0] = a0 + a1; stack[e0 + 1] = b0 + b1;
                        top = cbase + 1;
                    } else {
                        long sumA = 0, sumB = 0;
                        for (int j = 0; j < k; j++) {
                            int e = (cbase + j) * 2;
                            sumA += stack[e]; sumB += stack[e + 1];
                        }
                        for (int j = 0; j < k; j++) {
                            int e = (cbase + j) * 2;
                            long aj = stack[e], bj = stack[e + 1];
                            if (aj >= 2 && sumB != bj)
                                twoScore = twoScore.add(Int128.mulLong(aj * (aj - 1L), sumB - bj));
                            if (bj >= 2 && sumA != aj)
                                twoScore = twoScore.add(Int128.mulLong(bj * (bj - 1L), sumA - aj));
                        }
                        stack[cbase * 2] = (int) sumA; stack[cbase * 2 + 1] = (int) sumB;
                        top = cbase + 1;
                    }
                }
            }
        }
        return twoScore.halve();
    }

    /** CPU tree-walk path: build resident data once, then walk per split in parallel. */
    private void computeScoresCPUTreeWalk(List<BipartitionSplit> splitList, ClusterTable clusterTable,
                                          List<Tree> clusterTrees, List<Tree> partTrees,
                                          long[] scoreArray, double[] scoreArrayD, Int128[] scoreArrayI) {
        TreeWalkFrontier frontier = measureTreeWalkFrontier(partTrees);
        TreeWalkData d = buildTreeWalkData(splitList, clusterTable, clusterTrees, partTrees, frontier);
        int numSplits = splitList.size();
        int totalN = n;
        int cap = d.maxFrontier * 2;
        ThreadLocal<int[]> stackTL = ThreadLocal.withInitial(() -> new int[cap]);
        Logging.info("Weight table: CPU path (simple-tree-walk)  splits=%d  W=%d  clusters=%d  "
            + "trees=%d  tokens=%d  maxFrontier=%d (tree=%d, leaves=%d)",
            numSplits, d.W, d.numClusters, d.numTrees, d.nodeStream.length,
            d.maxFrontier, d.frontierTree, d.frontierLeaves);
        java.util.concurrent.atomic.AtomicInteger wDone = new java.util.concurrent.atomic.AtomicInteger(0);
        ProgressBar wBar = new ProgressBar("Scoring splits (CPU tree-walk)", numSplits);
        Threading.processRangeParallel(numSplits, idx -> {
            int[] stack = stackTL.get();
            if (useInt128)      scoreArrayI[idx] = computeScoreTreeWalkI(idx, d, totalN, stack);
            else if (useDouble) scoreArrayD[idx] = computeScoreTreeWalkD(idx, d, totalN, stack);
            else                scoreArray[idx]  = computeScoreTreeWalk(idx, d, totalN, stack);
            wBar.update(wDone.incrementAndGet());
        });
        wBar.done();
    }

    /** GPU tree-walk path: build resident data once, resolve batch size, call the native kernel. */
    private boolean computeScoresGPUTreeWalk(List<BipartitionSplit> splitList, ClusterTable clusterTable,
                                             List<Tree> clusterTrees, List<Tree> partTrees,
                                             long[] scoreArray, double[] scoreArrayD, Int128[] scoreArrayI) {
        // Early feasibility check against the ACTUAL postorder evaluation frontier,
        // not total taxa. A balanced tree with thousands of leaves can require only
        // a few dozen stack entries, while a pathologically ordered caterpillar can
        // still approach its leaf count.
        TreeWalkFrontier frontier = measureTreeWalkFrontier(partTrees);
        int selectedStackCap = frontier.entries <= 32 ? 32
            : frontier.entries <= 64 ? 64
            : frontier.entries <= 128 ? 128
            : frontier.entries <= 256 ? 256 : TW_GPU_STACK_CAP;
        Logging.info("GPU tree-walk frontier: required=%d entries (%d B/thread logically; "
            + "treeIndex=%d, leaves=%d), selected cap=%d (%d B/thread; compiled max=%d), numTaxa=%d",
            frontier.entries, frontier.entries * 2 * Integer.BYTES,
            frontier.treeIndex, frontier.leafCount, selectedStackCap,
            selectedStackCap * 2 * Integer.BYTES, TW_GPU_STACK_CAP, n);
        if (frontier.entries > TW_GPU_STACK_CAP) {
            Logging.info("GPU tree-walk: measured frontier=%d exceeds stack cap %d — using CPU tree-walk",
                frontier.entries, TW_GPU_STACK_CAP);
            return false;
        }

        Config cfg = Config.getInstance();
        TreeWalkData d = buildTreeWalkData(splitList, clusterTable, clusterTrees, partTrees, frontier);
        int numSplits = splitList.size();

        long residentMem = 8L * d.clusterBits.length
                         + 4L * ((long) d.nodeStream.length + d.treeNodeOffset.length);
        long perSplit = 4L * Integer.BYTES + 2L * d.W * Long.BYTES
                      + (useInt128 ? 2L : 1L) * Long.BYTES;

        int batchSizeHint; String batchDesc;
        if (!cfg.isGpuBatch()) {
            batchSizeHint = -1; batchDesc = "off (single launch)";
        } else if (cfg.getGpuNumBatches() > 0) {
            int N = cfg.getGpuNumBatches();
            batchSizeHint = (numSplits + N - 1) / N;
            batchDesc = N + " batches → batchSize=" + batchSizeHint;
        } else if (cfg.getGpuBatchSize() > 0) {
            batchSizeHint = cfg.getGpuBatchSize();
            batchDesc = "explicit batchSize=" + batchSizeHint;
        } else if (cfg.isGpuVramControlFactorSet()) {
            double F = cfg.getGpuVramControlFactor();
            long batchMem = (long) (F * residentMem);
            batchSizeHint = (int) Math.max(1, Math.min(numSplits, batchMem / perSplit));
            batchDesc = String.format("vram-control-factor=%.3f  resident=%.1f MB  batch→%d",
                F, residentMem / 1e6, (numSplits + batchSizeHint - 1) / Math.max(1, batchSizeHint));
        } else {
            long capBytes = (long) cfg.getGpuTreeWalkVramCapMiB() * 1024L * 1024L;
            batchSizeHint = (int) Math.max(1, Math.min((long) numSplits, capBytes / perSplit));
            int numBatches = (numSplits + batchSizeHint - 1) / batchSizeHint;
            batchDesc = String.format("auto scratch cap=%d MiB → %d batch%s (batchSize=%d)",
                cfg.getGpuTreeWalkVramCapMiB(), numBatches, numBatches == 1 ? "" : "es", batchSizeHint);
            if (numBatches > 1) {
                Logging.info("GPU tree-walk scratch capped at %d MiB: %d batches; "
                    + "raise --gpu-treewalk-vram-cap-mb only if this weight step is unusually launch-overhead limited",
                    cfg.getGpuTreeWalkVramCapMiB(), numBatches);
            }
        }

        Logging.info("Weight table: GPU path (simple-tree-walk)  splits=%d  W=%d  clusters=%d  trees=%d  tokens=%d  resident=%.1f MB  batching=%s",
            numSplits, d.W, d.numClusters, d.numTrees, d.nodeStream.length, residentMem / 1e6, batchDesc);

        long t1 = System.nanoTime();
        long[] twoScores = GPUWeightCalculator.computeWeightsTreeWalkGPU(
            d.splitCid, d.clusterBits, d.nodeStream, d.treeNodeOffset,
            numSplits, d.numClusters, d.numTrees, d.W, n, d.maxFrontier,
            batchSizeHint, cfg.getGpuVramFraction(), nativeScoreMode(),
            cfg.getGpuProgressIntervalSec());
        long gpuMs = (System.nanoTime() - t1) / 1_000_000;

        if (twoScores == null) {
            Logging.info("  GPU tree-walk kernel returned null after %d ms (infeasible)", gpuMs);
            return false;
        }
        Logging.info("  GPU tree-walk kernel returned in %d ms", gpuMs);
        unpackTwoScores(twoScores, scoreArray, scoreArrayD, scoreArrayI, numSplits);
        return true;
    }

    // -------------------------------------------------------------------------
    // CPU path (also the fallback when the GPU path is infeasible)
    // -------------------------------------------------------------------------

    private void computeScoresCPU(List<BipartitionSplit> splitList,
                                   Collection<PartitionTable.Entry> partitions,
                                   ClusterTable clusterTable,
                                   List<Tree> clusterTrees, List<Tree> partTrees,
                                   long[] scoreArray, double[] scoreArrayD, Int128[] scoreArrayI) {
        int numSplits = splitList.size();
        // CPU: parallel over splits (TRACE: single-threaded for deterministic output)
        if (Logging.isTrace()) {
            for (int idx = 0; idx < numSplits; idx++) {
                BipartitionSplit sp = splitList.get(idx);
                Logging.trace("SPLIT sz=%d|%d  lo=%s  hi=%s",
                    sp.lo.size, sp.hi.size, sp.lo, sp.hi);
                if (useInt128) {
                    scoreArrayI[idx] = computeScoreI(sp, partitions, clusterTable, clusterTrees, partTrees);
                    Logging.trace("  => score=%s", scoreArrayI[idx]);
                } else if (useDouble) {
                    scoreArrayD[idx] = computeScoreD(sp, partitions, clusterTable, clusterTrees, partTrees);
                    Logging.trace("  => score=%.6e", scoreArrayD[idx]);
                } else {
                    scoreArray[idx] = computeScore(sp, partitions, clusterTable, clusterTrees, partTrees);
                    Logging.trace("  => score=%d", scoreArray[idx]);
                }
            }
        } else {
            java.util.concurrent.atomic.AtomicInteger wDone = new java.util.concurrent.atomic.AtomicInteger(0);
            ProgressBar wBar = new ProgressBar("Scoring splits (CPU)", numSplits);
            Threading.processRangeParallel(numSplits, idx -> {
                if (useInt128) {
                    scoreArrayI[idx] = computeScoreI(splitList.get(idx), partitions, clusterTable,
                                                     clusterTrees, partTrees);
                } else if (useDouble) {
                    scoreArrayD[idx] = computeScoreD(splitList.get(idx), partitions, clusterTable,
                                                     clusterTrees, partTrees);
                } else {
                    scoreArray[idx] = computeScore(splitList.get(idx), partitions, clusterTable,
                                                   clusterTrees, partTrees);
                }
                wBar.update(wDone.incrementAndGet());
            });
            wBar.done();
        }
    }

    // -------------------------------------------------------------------------
    // Cluster-side intersection dispatch (single-range fast path / multi-range).
    //
    // For a single-range cluster (los == null) these are byte-identical to the
    // original IntersectionCounter calls.  For a multi-range cluster they sum the
    // intersection over the cluster's disjoint ranges (see multi-range-cluster-design.md
    // §5.1).  Centralizing here keeps all three numeric modes consistent.
    // -------------------------------------------------------------------------

    /** |M_range ∩ cluster| where cluster c (in tree tC) may be single- or multi-range. */
    private static int clusterIntersect(Tree tGT, int loGT, int hiGT,
                                        Tree tC, Cluster c, int sizeGTRange) {
        if (c.los != null)
            return IntersectionCounter.intersectMulti(tGT, loGT, hiGT, tC, c.los, c.his,
                                                      c.complement, sizeGTRange);
        return IntersectionCounter.intersect(tGT, loGT, hiGT, tC, c.left, c.right,
                                             c.complement, sizeGTRange);
    }

    private long computeScore(BipartitionSplit split,
                               Collection<PartitionTable.Entry> partitions,
                               ClusterTable clusterTable,
                               List<Tree> clusterTrees, List<Tree> partTrees) {
        // Retrieve exemplars for A (lo half) and B (hi half)
        ClusterTable.Entry eA = clusterTable.get(split.lo);
        ClusterTable.Entry eB = clusterTable.get(split.hi);
        if (eA == null || eB == null) return 0L;

        Cluster cA = eA.exemplar;
        Cluster cB = eB.exemplar;
        // Cluster positions are in completed trees; use clusterTrees for position lookup.
        Tree tA = clusterTrees.get(cA.treeIndex);
        Tree tB = clusterTrees.get(cB.treeIndex);
        int sizeA = cA.size;
        int sizeB = cB.size;
        int sizeC = n - sizeA - sizeB;
        if (sizeC < 0) return 0L;  // sanity

        long twoScore = 0L;

        for (PartitionTable.Entry pe : partitions) {
            Partition p = pe.exemplar;
            // Partition positions are in original trees; use partTrees for gene-tree lookup.
            Tree tGT = partTrees.get(p.treeIndex);

            // Polytomous rooted partition (d > 3): O(d) triplet formula.
            if (p.d > 3) {
                int[][] parts = polyChildIntersections(p, tGT, tA, cA, tB, cB);
                long twoQI = polyTwoQILong(parts[0], parts[1], p.d - 1);
                twoScore += (long) pe.frequency * twoQI;
                continue;
            }

            // M1 = [leftStart, leftEnd), M2 = [rightStart, rightEnd)
            int lo1 = p.leftStart,  hi1 = p.leftEnd;   // M1 range
            int lo2 = p.rightStart, hi2 = p.rightEnd;  // M2 range
            int sz1 = p.size1, sz2 = p.size2;

            // 4 core intersections
            int a0 = clusterIntersect(tGT, lo1, hi1, tA, cA, sz1);
            int a1 = clusterIntersect(tGT, lo2, hi2, tA, cA, sz2);
            int b0 = clusterIntersect(tGT, lo1, hi1, tB, cB, sz1);
            int b1 = clusterIntersect(tGT, lo2, hi2, tB, cB, sz2);

            long twoQI = computeTwoQI(a0, a1, b0, b1);
            Logging.trace("    CHILD-PART  tGT=%d sz=%d|%d a=[%d,%d] b=[%d,%d] "
                    + "2*tripletWeight=%d freq=%d",
                p.treeIndex, sz1, sz2, a0, a1, b0, b1, twoQI, pe.frequency);
            twoScore += (long) pe.frequency * twoQI;
        }

        // twoScore = 2 * score; every rooted-triplet weight is integral.
        return twoScore / 2;
    }

    /**
     * Compute twice the agreeing rooted-triplet count for child partitions
     * {@code A|B} and {@code M0|M1}.
     */
    private static long computeTwoQI(int a0, int a1, int b0, int b1) {
        long sameOrientation = (long) a0 * b1 * (a0 + b1 - 2L);
        long swappedOrientation = (long) a1 * b0 * (a1 + b0 - 2L);
        return sameOrientation + swappedOrientation;
    }

    // -------------------------------------------------------------------------
    // CPU path — DOUBLE variant (used when needsDoubleAccumulation() is true)
    //
    // Mirror of computeScore()/computeTwoQI() with floating-point accumulation so
    // very large taxon sets (where the exact integer 2·score overflows long) do
    // not wrap around.  The integer intersection matrix is computed identically;
    // only the triplet products and the running total are doubles.
    // -------------------------------------------------------------------------

    private double computeScoreD(BipartitionSplit split,
                                  Collection<PartitionTable.Entry> partitions,
                                  ClusterTable clusterTable,
                                  List<Tree> clusterTrees, List<Tree> partTrees) {
        ClusterTable.Entry eA = clusterTable.get(split.lo);
        ClusterTable.Entry eB = clusterTable.get(split.hi);
        if (eA == null || eB == null) return 0.0;

        Cluster cA = eA.exemplar;
        Cluster cB = eB.exemplar;
        Tree tA = clusterTrees.get(cA.treeIndex);
        Tree tB = clusterTrees.get(cB.treeIndex);
        int sizeA = cA.size;
        int sizeB = cB.size;
        int sizeC = n - sizeA - sizeB;
        if (sizeC < 0) return 0.0;

        double twoScore = 0.0;

        for (PartitionTable.Entry pe : partitions) {
            Partition p = pe.exemplar;
            Tree tGT = partTrees.get(p.treeIndex);

            if (p.d > 3) {   // polytomous: O(d) rooted weight (double accumulation)
                int[][] parts = polyChildIntersections(p, tGT, tA, cA, tB, cB);
                double twoQI = polyTwoQIDouble(parts[0], parts[1], p.d - 1);
                twoScore += (double) pe.frequency * twoQI;
                continue;
            }

            int lo1 = p.leftStart,  hi1 = p.leftEnd;
            int lo2 = p.rightStart, hi2 = p.rightEnd;
            int sz1 = p.size1, sz2 = p.size2;

            int a0 = clusterIntersect(tGT, lo1, hi1, tA, cA, sz1);
            int a1 = clusterIntersect(tGT, lo2, hi2, tA, cA, sz2);
            int b0 = clusterIntersect(tGT, lo1, hi1, tB, cB, sz1);
            int b1 = clusterIntersect(tGT, lo2, hi2, tB, cB, sz2);

            double twoQI = computeTwoQIDouble(a0, a1, b0, b1);
            twoScore += (double) pe.frequency * twoQI;
        }

        return twoScore / 2.0;
    }

    /** Floating-point mirror of {@link #computeTwoQI}; same formula, double accumulation. */
    private static double computeTwoQIDouble(int a0, int a1, int b0, int b1) {
        return (double) a0 * b1 * (a0 + b1 - 2.0)
             + (double) a1 * b0 * (a1 + b0 - 2.0);
    }

    // -------------------------------------------------------------------------
    // CPU path — INT128 variant (exact, overflow-free; used when the configured
    // large-score type is INT128).  Same intersection matrix as computeScore();
    // only the triplet products and the running total are 128-bit.
    // -------------------------------------------------------------------------

    private Int128 computeScoreI(BipartitionSplit split,
                                  Collection<PartitionTable.Entry> partitions,
                                  ClusterTable clusterTable,
                                  List<Tree> clusterTrees, List<Tree> partTrees) {
        ClusterTable.Entry eA = clusterTable.get(split.lo);
        ClusterTable.Entry eB = clusterTable.get(split.hi);
        if (eA == null || eB == null) return Int128.ZERO;

        Cluster cA = eA.exemplar;
        Cluster cB = eB.exemplar;
        Tree tA = clusterTrees.get(cA.treeIndex);
        Tree tB = clusterTrees.get(cB.treeIndex);
        int sizeA = cA.size;
        int sizeB = cB.size;
        int sizeC = n - sizeA - sizeB;
        if (sizeC < 0) return Int128.ZERO;

        Int128 twoScore = Int128.ZERO;

        for (PartitionTable.Entry pe : partitions) {
            Partition p = pe.exemplar;
            Tree tGT = partTrees.get(p.treeIndex);

            if (p.d > 3) {   // polytomous: O(d) rooted weight (exact 128-bit)
                int[][] parts = polyChildIntersections(p, tGT, tA, cA, tB, cB);
                Int128 twoQI = polyTwoQIInt128(parts[0], parts[1], p.d - 1);
                twoScore = twoScore.add(twoQI.mulScalar(pe.frequency));
                continue;
            }

            int lo1 = p.leftStart,  hi1 = p.leftEnd;
            int lo2 = p.rightStart, hi2 = p.rightEnd;
            int sz1 = p.size1, sz2 = p.size2;

            int a0 = clusterIntersect(tGT, lo1, hi1, tA, cA, sz1);
            int a1 = clusterIntersect(tGT, lo2, hi2, tA, cA, sz2);
            int b0 = clusterIntersect(tGT, lo1, hi1, tB, cB, sz1);
            int b1 = clusterIntersect(tGT, lo2, hi2, tB, cB, sz2);

            Int128 twoQI = computeTwoQIInt128(a0, a1, b0, b1);
            twoScore = twoScore.add(twoQI.mulScalar(pe.frequency));
        }

        return twoScore.halve();
    }

    /** Exact 128-bit mirror of {@link #computeTwoQI}. */
    private static Int128 computeTwoQIInt128(int a0, int a1, int b0, int b1) {
        Int128 same = (a0 == 0 || b1 == 0) ? Int128.ZERO
            : Int128.mulLong((long) a0 * b1, a0 + b1 - 2L);
        Int128 swapped = (a1 == 0 || b0 == 0) ? Int128.ZERO
            : Int128.mulLong((long) a1 * b0, a1 + b0 - 2L);
        return same.add(swapped);
    }

    // -------------------------------------------------------------------------
    // Rooted polytomy (d > 3), in O(d):
    //   2w = Σᵢ aᵢ(aᵢ-1)(Sb-bᵢ) + bᵢ(bᵢ-1)(Sa-aᵢ).
    // Only actual children participate; taxa outside the node are irrelevant.
    // -------------------------------------------------------------------------

    /**
     * Build the two intersection rows needed by the rooted-polytomy formula.
     */
    private static int[][] polyChildIntersections(Partition p, Tree tGT,
                                                   Tree tA, Cluster cA,
                                                   Tree tB, Cluster cB) {
        int childCount = p.d - 1;
        int[] a = new int[childCount], b = new int[childCount];
        for (int i = 0; i < childCount; i++) {
            int lo = p.partStarts[i], hi = p.partEnds[i], szi = p.sizes[i];
            a[i] = clusterIntersect(tGT, lo, hi, tA, cA, szi);
            b[i] = clusterIntersect(tGT, lo, hi, tB, cB, szi);
        }
        return new int[][]{ a, b };
    }

    private static long polyTwoQILong(int[] a, int[] b, int childCount) {
        long Sa = 0, Sb = 0;
        for (int i = 0; i < childCount; i++) { Sa += a[i]; Sb += b[i]; }
        long two = 0;
        for (int i = 0; i < childCount; i++) {
            long ai = a[i], bi = b[i];
            two += ai * (ai - 1) * (Sb - bi);
            two += bi * (bi - 1) * (Sa - ai);
        }
        return two;
    }

    private static double polyTwoQIDouble(int[] a, int[] b, int childCount) {
        double Sa = 0, Sb = 0;
        for (int i = 0; i < childCount; i++) { Sa += a[i]; Sb += b[i]; }
        double two = 0;
        for (int i = 0; i < childCount; i++) {
            double ai = a[i], bi = b[i];
            two += ai * (ai - 1) * (Sb - bi);
            two += bi * (bi - 1) * (Sa - ai);
        }
        return two;
    }

    private static Int128 polyTwoQIInt128(int[] a, int[] b, int childCount) {
        long Sa = 0, Sb = 0;
        for (int i = 0; i < childCount; i++) { Sa += a[i]; Sb += b[i]; }
        // Products can exceed signed long for extreme n, so form them in Int128.
        Int128 two = Int128.ZERO;
        for (int i = 0; i < childCount; i++) {
            long ai = a[i], bi = b[i];
            if (ai >= 2 && Sb != bi) two = two.add(Int128.mulLong(ai * (ai - 1), Sb - bi));
            if (bi >= 2 && Sa != ai) two = two.add(Int128.mulLong(bi * (bi - 1), Sa - ai));
        }
        return two;
    }

    // -------------------------------------------------------------------------
    // Queries
    // -------------------------------------------------------------------------

    /** The active score numeric type. */
    public Mode getMode()      { return mode; }
    public boolean isDouble()  { return useDouble; }
    public boolean isInt128()  { return useInt128; }

    /**
     * Score of a split as a long.  In DOUBLE/INT128 mode this returns the value
     * rounded/clamped to long (debug/verifier tooling only); the DP must use the
     * mode-matching accessor ({@link #getScoreD} / {@link #getScoreI}).
     */
    public long getScore(BipartitionSplit split) {
        int i = split.scoreIndex();
        if (i < 0 || i >= scoredSplits.length || scoredSplits[i] != split) return 0L;
        if (useInt128) return Math.round(scoresI[i].toDouble());
        if (useDouble) return Math.round(scoresD[i]);
        return scores[i];
    }

    /** Score of a split as a double (valid in all modes; approximate for INT128). */
    public double getScoreD(BipartitionSplit split) {
        int i = split.scoreIndex();
        if (i < 0 || i >= scoredSplits.length || scoredSplits[i] != split) return 0.0;
        if (useInt128) return scoresI[i].toDouble();
        if (useDouble) return scoresD[i];
        return (double) scores[i];
    }

    /** Score of a split as an exact Int128 (valid in all modes). */
    public Int128 getScoreI(BipartitionSplit split) {
        int i = split.scoreIndex();
        if (i < 0 || i >= scoredSplits.length || scoredSplits[i] != split) return Int128.ZERO;
        if (useInt128) return scoresI[i];
        if (useDouble) return Int128.ofLong(Math.round(scoresD[i]));
        return Int128.ofLong(scores[i]);
    }

    public long   getMaxScore()    { return useInt128 ? Math.round(getMaxScoreD())
                                          : useDouble ? Math.round(maxScoreD)   : maxScore; }
    public long   getTotalScore()  { return useInt128 ? Math.round(getTotalScoreD())
                                          : useDouble ? Math.round(totalScoreD) : totalScore; }
    public double getMaxScoreD()   { return useInt128 ? (maxScoreI == null ? 0.0 : maxScoreI.toDouble())
                                          : useDouble ? maxScoreD   : (double) maxScore; }
    public double getTotalScoreD() { return useInt128 ? totalScoreI.toDouble()
                                          : useDouble ? totalScoreD : (double) totalScore; }
    public int    size()           { return scoredSplits.length; }

    /** Iterate all (split, score) pairs (LONG mode only; empty in DOUBLE/INT128 mode). */
    public Set<Map.Entry<BipartitionSplit, Long>> entries() {
        if (mode != Mode.LONG) return Collections.emptySet();
        Set<Map.Entry<BipartitionSplit, Long>> out = new LinkedHashSet<>(
            Math.max(16, scoredSplits.length * 2));
        for (int i = 0; i < scoredSplits.length; i++) {
            out.add(new AbstractMap.SimpleImmutableEntry<>(scoredSplits[i], scores[i]));
        }
        return out;
    }
}
