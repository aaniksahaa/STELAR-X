package stelarx.greedy;

import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

/**
 * Sequential UPGMA on a small {@code g × g} similarity matrix, used by the
 * polytomy resolution phase.
 *
 * Why not reuse {@link stelarx.completion.UPGMAClusterer}?  That class uses
 * {@code Threading.processRangeParallel} internally.  If we call it from
 * inside an outer polytomy-pool worker (which is also submitted to the
 * shared {@code Threading} executor), the inner {@code processRangeParallel}
 * queues sub-tasks AND blocks on a latch — but every worker thread is busy
 * holding its polytomy task, so the sub-tasks never start: deadlock.
 *
 * For polytomy sizes ({@code g ≤ √(50+25n)}, typically ≤ 31 in practice),
 * the sequential UPGMA is also faster than the parallel one — internal
 * threading overhead dominates the actual work.
 *
 * Output format matches {@link stelarx.completion.UPGMAClusterer#build}:
 *   - {@link TreeNode#left}/{@link TreeNode#right} for the dendrogram structure
 *   - {@link Tree#postorderArray} listing the original cluster indices in
 *     left-to-right post-order
 *   - per-node {@code rangeStart}/{@code rangeEnd} stamps into postorderArray
 */
public final class MiniUPGMA {

    private MiniUPGMA() {}

    /**
     * Run UPGMA on the flat {@code n × n} similarity matrix and return the
     * dendrogram.  Time: O(n³); for n ≤ 31 this is ≤ 30k ops, trivial.
     */
    public static Tree build(double[] sim, int n, int treeIndex) {
        if (n <= 0) throw new IllegalArgumentException("MiniUPGMA: n <= 0");
        if (n == 1) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId = 0;
            leaf.rangeStart = 0;
            leaf.rangeEnd   = 1;
            return new Tree(treeIndex, leaf, new int[]{0}, new int[]{0}, 1, 1);
        }

        TreeNode[] clusterRoot = new TreeNode[n];
        double[]    weight     = new double[n];
        boolean[]   active     = new boolean[n];
        for (int i = 0; i < n; i++) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId   = i;
            clusterRoot[i] = leaf;
            weight[i]      = 1.0;
            active[i]      = true;
        }

        // Materialize a mutable copy so we can update rows in place
        double[][] mat = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) mat[i][j] = sim[i * n + j];
        }

        for (int iter = 0; iter < n - 1; iter++) {
            // Find max-similarity active pair
            int bestI = -1, bestJ = -1;
            double bestS = -Double.MAX_VALUE;
            for (int i = 0; i < n; i++) {
                if (!active[i]) continue;
                for (int j = i + 1; j < n; j++) {
                    if (!active[j]) continue;
                    if (mat[i][j] > bestS) { bestS = mat[i][j]; bestI = i; bestJ = j; }
                }
            }
            if (bestI < 0) break;          // disconnected — shouldn't happen on a full matrix

            // Merge J into I
            TreeNode newNode = new TreeNode();
            newNode.left  = clusterRoot[bestI];
            newNode.right = clusterRoot[bestJ];
            clusterRoot[bestI].parent = newNode;
            clusterRoot[bestJ].parent = newNode;
            clusterRoot[bestI] = newNode;

            double wI = weight[bestI], wJ = weight[bestJ];
            weight[bestI] = wI + wJ;
            active[bestJ] = false;

            // Weighted-average update of row I
            for (int k = 0; k < n; k++) {
                if (k == bestI || !active[k]) continue;
                double newIK = (mat[bestI][k] * wI + mat[bestJ][k] * wJ) / (wI + wJ);
                mat[bestI][k] = newIK;
                mat[k][bestI] = newIK;
            }
        }

        TreeNode root = null;
        for (int i = 0; i < n; i++) if (active[i]) { root = clusterRoot[i]; break; }
        if (root == null) root = clusterRoot[0];

        int[] postArr = new int[n];
        int[] posMap  = new int[n];
        int[] pos = {0};
        stampRanges(root, postArr, pos);
        for (int i = 0; i < n; i++) posMap[postArr[i]] = i;

        return new Tree(treeIndex, root, postArr, posMap, n, n);
    }

    /**
     * Exact average-linkage UPGMA via the <b>nearest-neighbour-chain</b> (RNN)
     * algorithm — O(n²) time / O(n²) space, no approximation.  Produces the same
     * dendrogram (cluster set) as {@link #build}: the merge ORDER differs, but the
     * Lance–Williams average-linkage update used here is byte-identical to
     * {@code build}'s, and average linkage is <i>reducible</i>, so the chain
     * algorithm provably yields the same UPGMA clustering.  Used only for large
     * polytomies (g &gt; 31) under {@code --stepb-process-large-polytomies}, where
     * {@code build}'s O(n³) closest-pair scan would dominate.
     *
     * Determinism: "nearest" = strictly-greatest similarity, ties broken by
     * smallest cluster index (matching {@code build}'s row-major scan direction);
     * the chain is seeded with the smallest active index.  Left/right child order
     * (lower original index = left) does not affect the emitted cluster set.
     */
    public static Tree buildFast(double[] sim, int n, int treeIndex) {
        if (n <= 0) throw new IllegalArgumentException("MiniUPGMA: n <= 0");
        if (n <= 2) return build(sim, n, treeIndex);   // trivial; reuse exact path

        TreeNode[] clusterRoot = new TreeNode[n];
        double[]   weight       = new double[n];
        boolean[]  active       = new boolean[n];
        for (int i = 0; i < n; i++) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId   = i;
            clusterRoot[i] = leaf;
            weight[i]      = 1.0;
            active[i]      = true;
        }

        // Operate on the caller's flat `sim` array IN PLACE (mat[i][j] == sim[i*n+j]).
        // This mutates `sim`; both call sites (Step A groupSim, Step B inducedSim)
        // discard it afterward.  Avoids a second d×d matrix (the old `double[n][n]`
        // copy) — bit-identical result, just no copy.
        double[] mat = sim;

        int[] chain   = new int[n + 1];
        int   chainSz = 0;
        int   remaining = n;
        int   nextSeed = 0;        // smallest-index seed scan cursor

        while (remaining > 1) {
            if (chainSz == 0) {
                while (!active[nextSeed]) nextSeed++;
                chain[chainSz++] = nextSeed;
            }
            int a = chain[chainSz - 1];
            int aRow = a * n;

            // nearest neighbour of a: max similarity, tie-break smallest index
            int b = -1;
            double bestS = -Double.MAX_VALUE;
            for (int c = 0; c < n; c++) {
                if (c == a || !active[c]) continue;
                if (mat[aRow + c] > bestS) { bestS = mat[aRow + c]; b = c; }
            }

            if (chainSz >= 2 && b == chain[chainSz - 2]) {
                // a and b are reciprocal nearest neighbours → merge
                chainSz -= 2;                       // pop a and b
                int lo = Math.min(a, b), hi = Math.max(a, b);
                double wLo = weight[lo], wHi = weight[hi];

                TreeNode newNode = new TreeNode();
                newNode.left  = clusterRoot[lo];
                newNode.right = clusterRoot[hi];
                clusterRoot[lo].parent = newNode;
                clusterRoot[hi].parent = newNode;
                clusterRoot[lo] = newNode;

                weight[lo] = wLo + wHi;
                active[hi] = false;
                int loRow = lo * n, hiRow = hi * n;
                for (int k = 0; k < n; k++) {
                    if (k == lo || !active[k]) continue;
                    double newLoK = (mat[loRow + k] * wLo + mat[hiRow + k] * wHi) / (wLo + wHi);
                    mat[loRow + k] = newLoK;
                    mat[k * n + lo] = newLoK;
                }
                remaining--;
            } else {
                chain[chainSz++] = b;               // extend the chain
            }
        }

        TreeNode root = null;
        for (int i = 0; i < n; i++) if (active[i]) { root = clusterRoot[i]; break; }
        if (root == null) root = clusterRoot[0];

        int[] postArr = new int[n];
        int[] posMap  = new int[n];
        int[] pos = {0};
        stampRanges(root, postArr, pos);
        for (int i = 0; i < n; i++) posMap[postArr[i]] = i;

        return new Tree(treeIndex, root, postArr, posMap, n, n);
    }

    private static void stampRanges(TreeNode node, int[] postArr, int[] pos) {
        int lo = pos[0];
        if (node.left == null) {           // leaf
            postArr[pos[0]++] = node.taxonId;
        } else {
            stampRanges(node.left,  postArr, pos);
            stampRanges(node.right, postArr, pos);
        }
        node.rangeStart = lo;
        node.rangeEnd   = pos[0];
    }
}
