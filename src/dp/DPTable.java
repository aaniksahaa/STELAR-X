package dp;

import utils.Logging;
import cluster.ClusterHash;
import cluster.ClusterTable;
import hash.PrefixHashArrays;
import tree.Tree;
import tree.TreeNode;

import java.util.*;

/**
 * DP search space: maps each cluster hash to its set of candidate bipartition splits.
 *
 * Built via tree-local transitions only, O(nk):
 *
 *   Type 1 -- for every internal node u (incl. root):
 *             sub(u) -> sub(left(u)) | sub(right(u))
 *
 *   Type 2 -- for every non-root internal node u whose parent is also non-root:
 *             [S \ sub(u)] -> sub(sibling(u)) | [S \ sub(parent(u))]
 *
 * For complete trees (Lg == S), Type 2 adds S\sub(root)=empty and is skipped.
 *
 * The root of the DP is the all-taxa cluster (from ClusterTable.getAllTaxaHash()).
 */
public class DPTable {

    // transitions[parentHash] = set of distinct BipartitionSplits for that cluster
    private final Map<ClusterHash, Set<BipartitionSplit>> transitions = new HashMap<>();

    private final ClusterHash rootHash; // allTaxaHash -- starting point of the DP
    private final int m;
    private final int n; // total taxa count

    private int totalEmitted = 0;
    private int uniqueSplits = 0;

    // -------------------------------------------------------------------------

    public DPTable(List<Tree> trees, PrefixHashArrays pref, ClusterTable clusterTable) {
        long t0 = System.nanoTime();
        this.m        = pref.numSeeds();
        this.rootHash = clusterTable.getAllTaxaHash();
        this.n        = rootHash.size;

        for (Tree tree : trees) {
            extractFromTree(tree, pref);
        }

        for (Set<BipartitionSplit> s : transitions.values()) uniqueSplits += s.size();

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("DP table: %d clusters with splits, %d unique splits (%d emitted) in %d ms",
            transitions.size(), uniqueSplits, totalEmitted, ms);
    }

    // -------------------------------------------------------------------------

    private void extractFromTree(Tree tree, PrefixHashArrays pref) {
        int ti = tree.treeIndex;
        int L  = tree.leafCount;

        // Annotate node ranges
        int nodeCount = tree.nodes != null ? tree.nodes.size() : 0;
        int[][] ranges = new int[Math.max(nodeCount, 1)][2];
        for (int[] r : ranges) { r[0] = -1; r[1] = -1; }
        annotateNode(tree.root, tree, ranges);

        emit(tree.root, ti, pref, ranges, /*parentRange=*/null);

        // Type 3: for incomplete trees add S -> Lg | (S\Lg)
        if (!tree.isComplete) {
            ClusterHash hLg  = hashRange(ti, 0, L, false, pref);
            ClusterHash hSLg = hashRange(ti, 0, L, true,  pref);
            if (hSLg.size > 0) {
                addTransition(rootHash, hLg, hSLg);
            }
        }
    }

    /**
     * Post-order recursion: emit Type 1 and Type 2 transitions for each internal node.
     */
    private void emit(TreeNode u, int ti, PrefixHashArrays pref, int[][] ranges,
                      int[] parentRange) {
        if (u.isLeaf()) return;
        if (u.childs == null || u.childs.size() != 2) return;

        TreeNode left  = u.childs.get(0);
        TreeNode right = u.childs.get(1);

        // Get this node's range
        int[] uRange = (u.index >= 0 && u.index < ranges.length) ? ranges[u.index] : null;
        if (uRange == null || uRange[0] < 0) return;

        int[] leftRange  = (left.index >= 0  && left.index  < ranges.length) ? ranges[left.index]  : null;
        int[] rightRange = (right.index >= 0 && right.index < ranges.length) ? ranges[right.index] : null;

        // Recurse first (post-order)
        emit(left,  ti, pref, ranges, uRange);
        emit(right, ti, pref, ranges, uRange);

        if (leftRange == null || leftRange[0] < 0) return;
        if (rightRange == null || rightRange[0] < 0) return;

        // Type 1: sub(u) -> sub(left) | sub(right)
        ClusterHash hU     = hashRange(ti, uRange[0],     uRange[1],     false, pref);
        ClusterHash hLeft  = hashRange(ti, leftRange[0],  leftRange[1],  false, pref);
        ClusterHash hRight = hashRange(ti, rightRange[0], rightRange[1], false, pref);
        addTransition(hU, hLeft, hRight);

        // Type 2: S\sub(u) -> sub(sibling) | S\sub(parent)
        // Only if parent exists (u is not root) -- handled by parentRange != null
        if (parentRange != null && parentRange[0] >= 0) {
            ClusterHash hCompParent = hashRange(ti, parentRange[0], parentRange[1], true, pref);
            if (hCompParent.size > 0) {
                // Type 2 for left child: S\sub(left) -> sub(right) | S\sub(u)
                ClusterHash hCompLeft = hashRange(ti, leftRange[0],  leftRange[1],  true, pref);
                addTransition(hCompLeft, hRight, hCompParent);

                // Type 2 for right child: S\sub(right) -> sub(left) | S\sub(u)
                ClusterHash hCompRight = hashRange(ti, rightRange[0], rightRange[1], true, pref);
                addTransition(hCompRight, hLeft, hCompParent);
            }
        }
    }

    // -------------------------------------------------------------------------

    private void addTransition(ClusterHash parent, ClusterHash a, ClusterHash b) {
        totalEmitted++;
        BipartitionSplit split = new BipartitionSplit(a, b);
        transitions.computeIfAbsent(parent, k -> new LinkedHashSet<>()).add(split);
    }

    private ClusterHash hashRange(int ti, int lo, int hi, boolean complement,
                                   PrefixHashArrays pref) {
        long[] rawSums = new long[m], rawXors = new long[m];
        for (int s = 0; s < m; s++) {
            rawSums[s] = complement ? pref.superCompSum(ti, s, lo, hi) : pref.rangeSum(ti, s, lo, hi);
            rawXors[s] = complement ? pref.superCompXor(ti, s, lo, hi) : pref.rangeXor(ti, s, lo, hi);
        }
        int sz = complement ? (n - (hi - lo)) : (hi - lo);
        return new ClusterHash(rawSums, rawXors, sz, m);
    }

    private int[] annotateNode(TreeNode node, Tree tree, int[][] ranges) {
        if (node == null || node.index < 0 || node.index >= ranges.length) return null;

        if (node.isLeaf()) {
            if (node.taxon != null) {
                int pos = tree.positionMap[node.taxon.id];
                if (pos >= 0) {
                    ranges[node.index][0] = pos;
                    ranges[node.index][1] = pos + 1;
                    return ranges[node.index];
                }
            }
            return null;
        }

        int lo = Integer.MAX_VALUE, hi = Integer.MIN_VALUE;
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                int[] cr = annotateNode(child, tree, ranges);
                if (cr != null && cr[0] >= 0) {
                    lo = Math.min(lo, cr[0]);
                    hi = Math.max(hi, cr[1]);
                }
            }
        }

        if (lo == Integer.MAX_VALUE) return null;
        ranges[node.index][0] = lo;
        ranges[node.index][1] = hi;
        return ranges[node.index];
    }

    // -------------------------------------------------------------------------
    // Queries
    // -------------------------------------------------------------------------

    public ClusterHash getRootHash()                       { return rootHash; }
    public Set<BipartitionSplit> getSplits(ClusterHash h)  { return transitions.getOrDefault(h, Collections.emptySet()); }
    public boolean hasSplits(ClusterHash h)                { return transitions.containsKey(h); }
    public int numClusters()                               { return transitions.size(); }
    public int numUniqueSplits()                           { return uniqueSplits; }
    public int numEmitted()                                { return totalEmitted; }
    public Set<Map.Entry<ClusterHash, Set<BipartitionSplit>>> entries() { return transitions.entrySet(); }
}
