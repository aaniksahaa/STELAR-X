package stelarx.completion;

import stelarx.Config;
import stelarx.Logging;
import stelarx.gpu.GPUDistanceMatrix;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.ProgressBar;
import stelarx.util.Threading;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Builds the average pairwise topological distance matrix from gene trees.
 *
 * CPU path  (buildCPU):
 *   For each gene tree T, for each internal node u at depth d_u:
 *     for every pair (a ∈ left(u), b ∈ right(u)):
 *       distSum[a][b] += depth(a) + depth(b) - 2*d_u
 *       cooccurrence[a][b]++
 *   O(k × n²) time, O(n²) space.
 *
 * GPU path  (buildGPU):
 *   1. Build Euler tour + sparse-table RMQ for every tree on CPU (parallel, O(kn log n)).
 *   2. Flatten into contiguous short[] arrays padded to E_max.
 *   3. Call native CUDA kernel: Δ-tree batching × B×B pair tiling.
 *      GPU VRAM = O(Δ n log n + B²);  CPU RAM = O(n²).
 *   4. Normalize.
 */
public class DistanceMatrixBuilder {

    // ── Public entry points ───────────────────────────────────────────────────

    public static DistanceMatrix buildCPU(List<Tree> trees, int n) {
        DistanceMatrix dm = new DistanceMatrix(n);

        ProgressBar bar  = new ProgressBar("Building distance matrix (CPU)", trees.size());
        int done = 0;
        for (Tree tree : trees) {
            accumulateCPU(tree, dm);
            bar.update(++done);
        }
        bar.done();

        dm.normalize();
        return dm;
    }

    public static DistanceMatrix buildGPU(List<Tree> trees, int n) {
        int k = trees.size();
        Logging.info("Building Euler tour + sparse tables for %d trees (parallel CPU)", k);

        // ── Step 1: Build per-tree Euler tours in parallel ────────────────────
        EulerTourBuilder.TourData[] tours = new EulerTourBuilder.TourData[k];
        ProgressBar eulerBar = new ProgressBar("Euler tour + RMQ build", k);
        AtomicInteger eulerDone = new AtomicInteger(0);

        Threading.processRangeParallel(k, i -> {
            tours[i] = EulerTourBuilder.build(trees.get(i), n);
            eulerBar.update(eulerDone.incrementAndGet());
        });
        eulerBar.done();

        // ── Step 2: Compute flat-array layout constants ───────────────────────
        int eMaxRaw = 0, logMaxRaw = 0;
        for (EulerTourBuilder.TourData td : tours) {
            if (td.tourLen > eMaxRaw)  eMaxRaw  = td.tourLen;
            if (td.log     > logMaxRaw) logMaxRaw = td.log;
        }
        // Pad E_max to next power-of-2 so RMQ index arithmetic stays in-bounds
        int ePadded = 1;
        while (ePadded < eMaxRaw) ePadded <<= 1;
        final int E_max   = ePadded;
        final int LOG_max = logMaxRaw;

        Logging.info("  E_max=%d  LOG=%d  flat euler %.1f MB  flat sparse %.1f MB",
            E_max, LOG_max,
            (double)k * E_max * 2 / 1e6,
            (double)k * LOG_max * E_max * 2 / 1e6);

        // ── Step 3: Flatten into contiguous Java arrays ───────────────────────
        // eulerDepths [k × E_max]            short (int16)
        // sparseMin   [k × LOG_max × E_max]  short (int16)
        // leafDepth   [k × n]                short (int16)  -1 if absent
        // firstOcc    [k × n]                int             -1 if absent
        // eulerLen    [k]                    int
        long edSize  = (long)k * E_max;
        long spSize  = (long)k * LOG_max * E_max;
        long ldSize  = (long)k * n;
        long foSize  = (long)k * n;

        short[] eulerDepths = new short[(int)edSize];
        short[] sparseMin   = new short[(int)spSize];
        short[] leafDepth   = new short[(int)ldSize];
        int[]   firstOcc    = new int  [(int)foSize];
        int[]   eulerLen    = new int  [k];

        Arrays.fill(leafDepth, (short)-1);
        Arrays.fill(firstOcc,  -1);

        ProgressBar flatBar  = new ProgressBar("Flattening tour data", k);
        AtomicInteger flatDone = new AtomicInteger(0);

        Threading.processRangeParallel(k, i -> {
            EulerTourBuilder.TourData td = tours[i];
            int len = td.tourLen;

            // euler depths
            long edOff = (long)i * E_max;
            for (int p = 0; p < len; p++) eulerDepths[(int)(edOff + p)] = td.depths[p];
            // (remaining E_max - len entries stay 0 — never accessed for valid l,r)

            // sparse table (only levels 0..td.log-1 are valid)
            long spOff = (long)i * LOG_max * E_max;
            for (int lvl = 0; lvl < td.log; lvl++) {
                int rowLen = Math.max(0, len - (1 << lvl) + 1);
                long dst = spOff + (long)lvl * E_max;
                for (int p = 0; p < rowLen; p++) sparseMin[(int)(dst + p)] = td.sparseMin[lvl][p];
            }

            // leaf depth (taxon → depth in this tree)
            long ldOff = (long)i * n;
            for (int a = 0; a < n; a++) {
                int fo = td.firstOcc[a];
                if (fo >= 0) {
                    firstOcc   [(int)(ldOff + a)] = fo;
                    // depth at firstOcc position = eulerDepths at that position
                    leafDepth  [(int)(ldOff + a)] = td.depths[fo];
                }
            }

            eulerLen[i] = len;
            flatBar.update(flatDone.incrementAndGet());
        });
        flatBar.done();

        // ── Step 4: Call GPU kernel ───────────────────────────────────────────
        Config cfg = Config.getInstance();
        int tileSizeB        = cfg.getGpuDistTileSizeB();
        double progressInterval = cfg.getGpuDpProgressInterval();   // reuse same setting
        int progressMaxSteps    = cfg.getGpuDpProgressMaxSteps();

        DistanceMatrix dm = new DistanceMatrix(n);
        GPUDistanceMatrix.computeDistancesGPU(
            eulerDepths, sparseMin, leafDepth, firstOcc, eulerLen,
            k, n, E_max, LOG_max,
            tileSizeB, progressInterval, progressMaxSteps,
            dm.distSum, dm.cooccurrence
        );

        dm.normalize();
        return dm;
    }

    // ── CPU accumulation (per-tree) ───────────────────────────────────────────

    private static void accumulateCPU(Tree tree, DistanceMatrix dm) {
        int n = dm.n;
        int[] leafDepths = new int[n];
        Arrays.fill(leafDepths, -1);
        computeLeafDepths(tree.root, 0, leafDepths);
        accumulatePairs(tree.root, 0, leafDepths, tree, dm);
    }

    private static void computeLeafDepths(TreeNode node, int depth, int[] leafDepths) {
        if (node.isLeaf()) { leafDepths[node.taxonId] = depth; return; }
        computeLeafDepths(node.left,  depth + 1, leafDepths);
        computeLeafDepths(node.right, depth + 1, leafDepths);
    }

    private static void accumulatePairs(TreeNode node, int depth,
                                         int[] leafDepths, Tree tree, DistanceMatrix dm) {
        if (node.isLeaf()) return;
        accumulatePairs(node.left,  depth + 1, leafDepths, tree, dm);
        accumulatePairs(node.right, depth + 1, leafDepths, tree, dm);

        int n    = dm.n;
        int loL  = node.left.rangeStart,  hiL = node.left.rangeEnd;
        int loR  = node.right.rangeStart, hiR = node.right.rangeEnd;
        double twoDu = 2.0 * depth;

        for (int pi = loL; pi < hiL; pi++) {
            int    ta = tree.postorderArray[pi];
            double da = leafDepths[ta];
            for (int pj = loR; pj < hiR; pj++) {
                int    tb = tree.postorderArray[pj];
                double d  = da + leafDepths[tb] - twoDu;
                dm.distSum[ta * n + tb] += d;
                dm.distSum[tb * n + ta] += d;
                dm.cooccurrence[ta * n + tb]++;
                dm.cooccurrence[tb * n + ta]++;
            }
        }
    }
}
