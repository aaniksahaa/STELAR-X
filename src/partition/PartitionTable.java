package partition;

import utils.Logging;
import cluster.ClusterHash;
import hash.PrefixHashArrays;
import tree.Tree;
import tree.TreeNode;

import java.util.*;

/**
 * Table of unique gene-tree tripartitions with their frequencies.
 *
 * For each non-root internal node u of each rooted gene tree g we extract:
 *   part1 = sub(left(u))    range [L.start, L.end)   -- left subtree
 *   part2 = sub(right(u))   range [R.start, R.end)   -- right subtree
 *   part3 = Lg \ sub(u)     complement of [u.start, u.end) w.r.t. Lg
 *
 * The root node's two children define a bipartition with part3 = empty set,
 * which contributes 0 weight in triplet scoring -- so we skip the root.
 *
 * Deduplication: PartitionHash is order-invariant over (part1, part2).
 */
public class PartitionTable {

    public static final class Entry {
        public final PartitionHash  hash;
        public final RangePartition exemplar;
        public int                  frequency;

        Entry(PartitionHash h, RangePartition p) { this.hash = h; this.exemplar = p; this.frequency = 1; }
    }

    private final Map<PartitionHash, Entry> table = new HashMap<>();
    private final int m;

    // -------------------------------------------------------------------------

    public PartitionTable(List<Tree> trees, PrefixHashArrays pref) {
        long t0 = System.nanoTime();
        this.m = pref.numSeeds();

        int totalCandidates = 0;
        for (Tree tree : trees) {
            totalCandidates += extractFromTree(tree, pref);
        }

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Partition extraction: %d candidates -> %d unique tripartitions in %d ms",
            totalCandidates, table.size(), ms);
    }

    private int extractFromTree(Tree tree, PrefixHashArrays pref) {
        int ti = tree.treeIndex;
        int L  = tree.leafCount;
        int[] count = {0};

        // Annotate node ranges
        int nodeCount = tree.nodes != null ? tree.nodes.size() : 0;
        int[][] ranges = new int[Math.max(nodeCount, 1)][2];
        for (int[] r : ranges) { r[0] = -1; r[1] = -1; }
        annotateNode(tree.root, tree, ranges);

        extractNode(tree.root, ti, L, pref, count, ranges, /*isRoot=*/true);
        return count[0];
    }

    /**
     * Recurse post-order. For each non-root internal node u with exactly 2 children,
     * register the tripartition (left | right | complement_of_node_range_w.r.t_Lg).
     */
    private void extractNode(TreeNode node, int ti, int L,
                              PrefixHashArrays pref, int[] count,
                              int[][] ranges, boolean isRoot) {
        if (node.isLeaf()) return;

        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                extractNode(child, ti, L, pref, count, ranges, /*isRoot=*/false);
            }
        }

        if (isRoot) return;  // root gives bipartition with empty part3 -- skip

        // Need exactly 2 children for a binary tripartition
        if (node.childs == null || node.childs.size() != 2) return;

        TreeNode leftChild  = node.childs.get(0);
        TreeNode rightChild = node.childs.get(1);

        if (leftChild.index < 0 || leftChild.index >= ranges.length) return;
        if (rightChild.index < 0 || rightChild.index >= ranges.length) return;
        if (node.index < 0 || node.index >= ranges.length) return;

        int lStart = ranges[leftChild.index][0],  lEnd = ranges[leftChild.index][1];
        int rStart = ranges[rightChild.index][0], rEnd = ranges[rightChild.index][1];
        int pStart = ranges[node.index][0],       pEnd = ranges[node.index][1];

        if (lStart < 0 || lEnd <= lStart) return;
        if (rStart < 0 || rEnd <= rStart) return;
        if (pStart < 0 || pEnd <= pStart) return;

        int sz1 = lEnd - lStart;
        int sz2 = rEnd - rStart;
        int sz3 = L - (pEnd - pStart);   // complement size w.r.t. Lg

        if (sz3 == 0) return;  // only happens when node is root, already skipped

        // Build ClusterHash for each part
        ClusterHash h1 = buildHash(ti, lStart, lEnd, false, sz1, pref);
        ClusterHash h2 = buildHash(ti, rStart, rEnd, false, sz2, pref);
        ClusterHash h3 = buildHash(ti, pStart, pEnd, true,  sz3, pref); // complement w.r.t. Lg

        PartitionHash ph = new PartitionHash(h1, h2, h3);

        Entry existing = table.get(ph);
        if (existing != null) {
            existing.frequency++;
        } else {
            RangePartition p = new RangePartition(h1, h2, h3, sz1, sz2, sz3,
                                                  ti, lStart, lEnd, rStart, rEnd);
            table.put(ph, new Entry(ph, p));
        }
        count[0]++;
    }

    /** Compute a ClusterHash for a subtree range or its complement w.r.t. Lg. */
    private ClusterHash buildHash(int ti, int lo, int hi, boolean complement,
                                   int size, PrefixHashArrays pref) {
        long[] rawSums = new long[m], rawXors = new long[m];
        for (int s = 0; s < m; s++) {
            rawSums[s] = complement ? pref.compSum(ti, s, lo, hi)
                                    : pref.rangeSum(ti, s, lo, hi);
            rawXors[s] = complement ? pref.compXor(ti, s, lo, hi)
                                    : pref.rangeXor(ti, s, lo, hi);
        }
        return new ClusterHash(rawSums, rawXors, size, m);
    }

    /** Recursively annotate each node with its [rangeStart, rangeEnd) in postorderArray. */
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

    public Entry get(PartitionHash ph) { return table.get(ph); }
    public int size()                  { return table.size(); }
    public Collection<Entry> entries() { return table.values(); }
}
