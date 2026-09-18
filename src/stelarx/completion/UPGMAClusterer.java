package stelarx.completion;

import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.Threading;

import java.util.Arrays;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;

/**
 * Efficient UPGMA clusterer using flat primitive arrays + pointer-based tree construction.
 *
 * <h3>Algorithm summary</h3>
 * Standard UPGMA agglomerative clustering (n−1 merges):
 * <ol>
 *   <li>Find the active pair (I, J) with highest similarity.</li>
 *   <li>Merge J into I: build a new TreeNode parent, update weights.</li>
 *   <li>Update sim(I, k) for every other active k via weighted average.</li>
 *   <li>Repair per-row best-neighbor caches; recompute stale rows.</li>
 * </ol>
 * After n−1 merges the single remaining cluster root is the UPGMA dendrogram.
 * A post-order walk assigns rangeStart/rangeEnd/postorderArray and returns a Tree.
 *
 * <h3>Data structures</h3>
 * <ul>
 *   <li>{@code double[] mat} — upper triangle (i ≤ j), size n(n+1)/2.</li>
 *   <li>{@code double[] bestSim}, {@code int[] bestJ} — per-row best-neighbor cache,
 *       O(n) total, fits in L2 cache for n ≤ 25,000.</li>
 *   <li>{@code int[] activeArr} + {@code int[] posInActive} — swap-and-shrink active set:
 *       removing a deactivated cluster J from activeArr takes O(1) (swap with last
 *       element, decrement activeCount).  Eliminates the O(n²) rebuild overhead of
 *       a naïve "scan-all-n-entries-each-iteration" approach.</li>
 *   <li>{@code TreeNode[] clusterRoot} — pointer to current cluster's subtree root.
 *       No BitSets; bipartitions are extracted after the fact via rangeStart/rangeEnd.</li>
 * </ul>
 *
 * <h3>Parallelization (intra-iteration)</h3>
 * <ul>
 *   <li><b>Find-best</b>: parallel reduction over {@code bestSim[]} — each thread scans
 *       its chunk, keeps a local best, then one serial reduce over T local winners.</li>
 *   <li><b>Row update</b>: parallel-for over active k ≠ I.  Each k is owned by exactly
 *       one thread; writes to {@code mat[idx(I,k)]}, {@code bestSim[k]}, {@code bestJ[k]}
 *       are disjoint across threads — no synchronization needed on those fields.</li>
 *   <li><b>Stale-row recompute</b>: stale k indices collected in a
 *       {@code ConcurrentLinkedQueue} during the parallel update, then dispatched as a
 *       second parallel batch.</li>
 *   <li><b>Recompute row I</b>: always serial (I's entire row changed).</li>
 *   <li>The outer merge sequence is inherently serial; parallelism is intra-iteration.</li>
 * </ul>
 *
 * <h3>Complexity</h3>
 * <ul>
 *   <li>Time: O(n²) expected (O(1) cache updates per active k in the common case);
 *       O(n³) worst case (every row stale every iteration — degenerate input).</li>
 *   <li>Memory: O(n²) for the upper triangle; O(n) for all other structures.</li>
 * </ul>
 */
public class UPGMAClusterer {

    // ── Instance state ────────────────────────────────────────────────────────

    private final int      n;
    private final double[] mat;         // upper triangle, size n*(n+1)/2
    private final double[] bestSim;     // per-row best similarity
    private final int[]    bestJ;       // per-row best neighbor index
    private final int[]    weight;
    private final TreeNode[] clusterRoot;  // current subtree root per cluster

    // Active-set maintained via swap-and-shrink (O(1) removal)
    private final int[] activeArr;      // activeArr[0..activeCount-1] = active indices
    private final int[] posInActive;    // posInActive[i] = position of i in activeArr
    private int         activeCount;

    // ── Construction ─────────────────────────────────────────────────────────

    private UPGMAClusterer(double[] sim, int n) {
        this.n           = n;
        int triSize      = n * (n + 1) / 2;
        this.mat         = new double[triSize];
        this.bestSim     = new double[n];
        this.bestJ       = new int[n];
        this.weight      = new int[n];
        this.clusterRoot = new TreeNode[n];
        this.activeArr   = new int[n];
        this.posInActive = new int[n];
        this.activeCount = n;

        // Copy upper triangle from flat n×n sim array
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                mat[idx(i, j)] = sim[i * n + j];
            }
        }

        // Initialise active set (identity ordering) and leaf TreeNodes
        for (int i = 0; i < n; i++) {
            activeArr[i]   = i;
            posInActive[i] = i;
            weight[i]      = 1;
            TreeNode leaf  = new TreeNode();
            leaf.taxonId   = i;
            clusterRoot[i] = leaf;
        }

        // Build per-row best-neighbor cache (single sequential pass per row)
        for (int i = 0; i < n; i++) {
            double bs = Double.NEGATIVE_INFINITY;
            int    bj = -1;
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                double s = mat[idx(i, j)];
                if (s > bs) { bs = s; bj = j; }
            }
            bestSim[i] = bs;
            bestJ[i]   = bj;
        }
    }

    // ── Public entry point ────────────────────────────────────────────────────

    /**
     * Run UPGMA on {@code sim[]} (flat n×n, row-major) and return the resulting
     * dendrogram as a {@link Tree} with all n taxa.
     *
     * @param sim       flat n×n similarity matrix (row-major, {@code sim[i*n+j]})
     * @param n         number of taxa
     * @param treeIndex tree index to embed in the returned Tree object
     */
    public static Tree build(double[] sim, int n, int treeIndex) {
        if (n == 1) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId   = 0;
            leaf.rangeStart = 0;
            leaf.rangeEnd   = 1;
            return new Tree(treeIndex, leaf,
                new int[]{0}, new int[]{0}, 1, 1);
        }
        UPGMAClusterer u = new UPGMAClusterer(sim, n);
        return u.run(treeIndex);
    }

    /** Dispatch to the unchanged dense implementation or the exact large-N implementation. */
    public static Tree build(SimilarityMatrix sim, int treeIndex) {
        if (!sim.isPacked()) return build(sim.sim, sim.n, treeIndex);
        if (sim.n == 1) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId = 0;
            leaf.rangeStart = 0;
            leaf.rangeEnd = 1;
            return new Tree(treeIndex, leaf, new int[]{0}, new int[]{0}, 1, 1);
        }
        return SegmentedUPGMAClusterer.build(sim, treeIndex);
    }

    // ── Core UPGMA loop ───────────────────────────────────────────────────────

    private Tree run(int treeIndex) {
        int T = Threading.getNumThreads();

        for (int iter = 0; iter < n - 1; iter++) {

            // ── Step 1: Find global best I (parallel reduction) ───────────────
            int I = findBestI(T);
            int J = bestJ[I];

            // ── Step 2: Build new internal TreeNode (merge J into I) ──────────
            TreeNode newNode = new TreeNode();
            newNode.left              = clusterRoot[I];
            newNode.right             = clusterRoot[J];
            clusterRoot[I].parent     = newNode;
            clusterRoot[J].parent     = newNode;
            clusterRoot[I]            = newNode;

            final double wI = weight[I];
            final double wJ = weight[J];
            weight[I] += weight[J];

            // ── Step 3: Deactivate J in O(1) via swap-and-shrink ─────────────
            removeActive(J);

            // Snapshot for parallel phase (activeCount and activeArr just updated)
            final int snapCount = activeCount;

            // ── Step 4: Parallel update of sim(I, k) and row-k caches ─────────
            // Each k in activeArr[0..snapCount-1] is processed by exactly one thread.
            // Writes to mat[idx(I,k)], bestSim[k], bestJ[k] are disjoint across threads.
            final int finalI = I;
            final int finalJ = J;
            ConcurrentLinkedQueue<Integer> staleQueue = new ConcurrentLinkedQueue<>();

            Threading.processRangeParallel(snapCount, ai -> {
                int k = activeArr[ai];
                if (k == finalI) return;

                double oldIK = mat[idx(finalI, k)];
                double oldJK = mat[idx(finalJ, k)];
                double newIK = (oldIK * wI + oldJK * wJ) / (wI + wJ);
                mat[idx(finalI, k)] = newIK;

                // Update row-k cache
                if (bestJ[k] == finalJ) {
                    // J gone: forced full recompute
                    staleQueue.add(k);
                } else if (newIK > bestSim[k]) {
                    // I is now the best neighbor of k
                    bestJ[k]   = finalI;
                    bestSim[k] = newIK;
                } else if (bestJ[k] == finalI && newIK < bestSim[k]) {
                    // cached best got worse; something else may now be better
                    staleQueue.add(k);
                }
                // else: bestJ[k] is some active m ≠ I,J; sim(k,m) unchanged; valid
            });

            // ── Step 5: Parallel recompute of stale rows ──────────────────────
            Integer[] staleArr = staleQueue.toArray(new Integer[0]);
            if (staleArr.length > 0) {
                Threading.processRangeParallel(staleArr.length, si ->
                    recomputeBest(staleArr[si])
                );
            }

            // ── Step 6: Recompute row I (serial — its entire row changed) ─────
            recomputeBest(I);
        }

        return buildTree(treeIndex);
    }

    // ── Parallel find-best ────────────────────────────────────────────────────

    /**
     * Parallel reduction: find active i with maximum bestSim[i].
     * T threads each scan a disjoint chunk; one serial reduce over T local winners.
     */
    private int findBestI(int T) {
        // Snap the current active slice for the parallel scan
        final int snap = activeCount;
        int actual     = Math.min(T, snap);
        int chunk      = Math.max(1, (snap + actual - 1) / actual);

        int[]    localI = new int   [actual];
        double[] localS = new double[actual];
        Arrays.fill(localI, -1);
        Arrays.fill(localS, Double.NEGATIVE_INFINITY);

        CountDownLatch latch = new CountDownLatch(actual);
        for (int t = 0; t < actual; t++) {
            int tid = t;
            int lo  = t * chunk;
            int hi  = Math.min(lo + chunk, snap);
            Threading.submit(() -> {
                int    bI = -1;
                double bS = Double.NEGATIVE_INFINITY;
                for (int ai = lo; ai < hi; ai++) {
                    int i = activeArr[ai];
                    if (bestSim[i] > bS) { bS = bestSim[i]; bI = i; }
                }
                localI[tid] = bI;
                localS[tid] = bS;
                latch.countDown();
            });
        }
        try { latch.await(); }
        catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException(e);
        }

        int    I    = -1;
        double best = Double.NEGATIVE_INFINITY;
        for (int t = 0; t < actual; t++) {
            if (localS[t] > best) { best = localS[t]; I = localI[t]; }
        }
        return I;
    }

    // ── Active-set management (O(1) removal via swap-and-shrink) ─────────────

    /**
     * Remove cluster j from the active set in O(1).
     * Swaps j with the last active element, then decrements activeCount.
     */
    private void removeActive(int j) {
        int pos  = posInActive[j];
        int last = activeArr[activeCount - 1];
        // Put 'last' into the slot that 'j' occupied
        activeArr[pos]    = last;
        posInActive[last] = pos;
        // j is now logically removed
        activeCount--;
    }

    // ── Row best-neighbor recompute ───────────────────────────────────────────

    /** Full linear scan of active row k to recompute bestSim[k] / bestJ[k]. */
    private void recomputeBest(int k) {
        double bs = Double.NEGATIVE_INFINITY;
        int    bj = -1;
        for (int ai = 0; ai < activeCount; ai++) {
            int m = activeArr[ai];
            if (m == k) continue;
            double s = mat[idx(k, m)];
            if (s > bs) { bs = s; bj = m; }
        }
        bestSim[k] = bs;
        bestJ[k]   = bj;
    }

    // ── Upper-triangle index ──────────────────────────────────────────────────

    /** Index into the upper-triangle array (i ≤ j).  Swaps if needed. */
    private int idx(int i, int j) {
        if (i > j) { int t = i; i = j; j = t; }
        return i * n - i * (i + 1) / 2 + j;
    }

    // ── Tree construction ─────────────────────────────────────────────────────

    /**
     * After all merges, the single remaining active cluster's root is the UPGMA
     * dendrogram root.  Assigns rangeStart/rangeEnd + postorderArray via DFS,
     * then wraps everything in a {@link Tree}.
     */
    private Tree buildTree(int treeIndex) {
        // The one remaining active cluster
        TreeNode root = clusterRoot[activeArr[0]];

        int[]   postArr = new int[n];
        int[]   cursor  = {0};
        assignRanges(root, postArr, cursor);

        int[]   posMap  = new int[n];
        Arrays.fill(posMap, -1);
        for (int pos = 0; pos < n; pos++) posMap[postArr[pos]] = pos;

        return new Tree(treeIndex, root, postArr, posMap, n, n);
    }

    /**
     * Post-order DFS: fills postorderArray left-to-right and sets
     * rangeStart/rangeEnd on every node.
     */
    private void assignRanges(TreeNode node, int[] postArr, int[] cursor) {
        if (node.isLeaf()) {
            node.rangeStart = cursor[0];
            postArr[cursor[0]++] = node.taxonId;
            node.rangeEnd   = cursor[0];
        } else {
            node.rangeStart = cursor[0];
            assignRanges(node.left,  postArr, cursor);
            assignRanges(node.right, postArr, cursor);
            node.rangeEnd   = cursor[0];
        }
    }
}
