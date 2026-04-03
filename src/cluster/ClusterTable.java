package cluster;

import utils.Logging;
import hash.PrefixHashArrays;
import tree.Tree;
import tree.TreeNode;

import java.util.*;

/**
 * The cluster set X: all unique clusters extracted from gene trees.
 *
 * For each rooted gene tree we walk every internal node u
 * (excluding root) and register two clusters:
 *   1. sub(u)      -- the subtree range [u.left, u.right)
 *   2. S \ sub(u)  -- super-complement w.r.t. ALL taxa (if size > 0)
 *
 * Also registers singleton leaf clusters as DP base cases.
 * Empty clusters (size 0) are excluded.
 *
 * Deduplication is done by ClusterHash. One exemplar Cluster is kept per unique hash.
 * The table also maintains size-binned lists for DP space construction.
 */
public class ClusterTable {

    /** Entry in the cluster hash table. */
    public static final class Entry {
        public final ClusterHash hash;
        public final Cluster     exemplar;  // any one cluster with this taxa set
        public int               frequency; // how many times this exact taxa set appeared

        Entry(ClusterHash h, Cluster c) { this.hash = h; this.exemplar = c; this.frequency = 1; }
    }

    // Main table: ClusterHash -> Entry
    private final Map<ClusterHash, Entry> table = new HashMap<>();

    // Size bins: size -> list of ClusterHash objects of that size
    private final Map<Integer, List<ClusterHash>> sizeBins = new HashMap<>();

    // Special: the all-taxa cluster hash (DP root)
    private ClusterHash allTaxaHash;

    private final int m; // number of hash seeds

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    public ClusterTable(List<Tree> trees, PrefixHashArrays pref, int numTaxa) {
        long t0 = System.nanoTime();
        this.m = pref.numSeeds();

        // Build all-taxa hash from prefix arrays
        long[] atSums = new long[m], atXors = new long[m];
        for (int s = 0; s < m; s++) {
            atSums[s] = pref.allTaxaSum(s);
            atXors[s] = pref.allTaxaXor(s);
        }
        allTaxaHash = new ClusterHash(atSums, atXors, numTaxa, m);

        int totalCandidates = 0;
        for (Tree tree : trees) {
            totalCandidates += extractFromTree(tree, pref, numTaxa);
        }

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Cluster extraction: %d candidates -> %d unique clusters in %d ms",
            totalCandidates, table.size(), ms);
    }

    /**
     * Extract clusters from one tree using node-range annotation.
     * Returns the number of candidates generated (before dedup).
     */
    private int extractFromTree(Tree tree, PrefixHashArrays pref, int numTaxa) {
        int ti = tree.treeIndex;
        int L  = tree.leafCount;
        int[] count = {0};

        // Annotate every node with its [rangeStart, rangeEnd) in postorderArray
        int nodeCount = tree.nodes != null ? tree.nodes.size() : 0;
        int[][] ranges = new int[Math.max(nodeCount, 1)][2];
        for (int[] r : ranges) { r[0] = -1; r[1] = -1; }
        annotateNode(tree.root, tree, ranges);

        // Walk annotated tree and register clusters
        walkAnnotated(tree.root, ti, L, pref, numTaxa, count, ranges, /*isRoot=*/true);
        return count[0];
    }

    /**
     * Post-order walk using pre-annotated ranges.
     * For every node except root: register sub(u) and S\sub(u).
     */
    private void walkAnnotated(TreeNode node, int ti, int L,
                               PrefixHashArrays pref, int numTaxa, int[] count,
                               int[][] ranges, boolean isRoot) {
        if (!node.isLeaf() && node.childs != null) {
            for (TreeNode child : node.childs) {
                walkAnnotated(child, ti, L, pref, numTaxa, count, ranges, /*isRoot=*/false);
            }
        }

        if (isRoot) return;  // skip root -- it is the all-taxa cluster

        if (node.index < 0 || node.index >= ranges.length) return;
        int lo = ranges[node.index][0];
        int hi = ranges[node.index][1];
        if (lo < 0 || hi <= lo) return;

        int rangeSize = hi - lo;

        // 1. Subtree cluster [lo, hi)
        registerCluster(ti, lo, hi, false, rangeSize, L, pref, numTaxa);
        count[0]++;

        // 2. Super-complement: S \ [lo, hi) (w.r.t. ALL taxa)
        int superCompSize = numTaxa - rangeSize;
        if (superCompSize > 0) {
            registerCluster(ti, lo, hi, true, superCompSize, numTaxa, pref, numTaxa);
            count[0]++;
        }
    }

    /**
     * Recursively annotate each node with [rangeStart, rangeEnd) in postorderArray.
     * Returns the range for this node, or null if the node has no leaves in the tree.
     */
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

    /**
     * Compute hash for a cluster and insert (or increment frequency) in the table.
     */
    private void registerCluster(int ti, int lo, int hi, boolean complement,
                                  int size, int leafCount,
                                  PrefixHashArrays pref, int numTaxa) {
        long[] rawSums = new long[m], rawXors = new long[m];
        for (int s = 0; s < m; s++) {
            if (!complement) {
                rawSums[s] = pref.rangeSum(ti, s, lo, hi);
                rawXors[s] = pref.rangeXor(ti, s, lo, hi);
            } else {
                rawSums[s] = pref.superCompSum(ti, s, lo, hi);
                rawXors[s] = pref.superCompXor(ti, s, lo, hi);
            }
        }

        ClusterHash hash = new ClusterHash(rawSums, rawXors, size, m);

        // Skip the all-taxa cluster (it's the DP root, not in X)
        if (hash.equals(allTaxaHash)) return;

        Entry existing = table.get(hash);
        if (existing != null) {
            existing.frequency++;
        } else {
            Cluster exemplar = new Cluster(ti, lo, hi, complement, leafCount);
            Entry entry = new Entry(hash, exemplar);
            table.put(hash, entry);
            sizeBins.computeIfAbsent(size, k -> new ArrayList<>()).add(hash);
        }
    }

    // -------------------------------------------------------------------------
    // Queries
    // -------------------------------------------------------------------------

    public Entry get(ClusterHash hash)      { return table.get(hash); }
    public boolean contains(ClusterHash h)  { return table.containsKey(h); }
    public int size()                       { return table.size(); }
    public ClusterHash getAllTaxaHash()     { return allTaxaHash; }
    public Collection<Entry> entries()     { return table.values(); }
    public int numSeeds()                  { return m; }

    /** All cluster hashes of a given size. */
    public List<ClusterHash> getBySize(int sz) {
        return sizeBins.getOrDefault(sz, Collections.emptyList());
    }

    /** All sizes present in X. */
    public Set<Integer> sizes() { return sizeBins.keySet(); }
}
