package stelarx.cluster;

import stelarx.Config;
import stelarx.Logging;
import stelarx.hash.PrefixHashArrays;
import stelarx.util.ProgressBar;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.*;

/**
 * The cluster set X: all unique clusters extracted from gene trees.
 *
 * For each rooted gene tree, we walk every node u (excluding the root) and
 * register its descendant cluster sub(u). Complements are not rooted clades.
 *
 * Also registers the all-taxa cluster (DP root) separately.
 * Singleton clusters (size 1) are included -- they are DP base cases.
 * Empty clusters (size 0) are excluded.
 *
 * Deduplication is done by ClusterHash. One exemplar Cluster is kept per unique hash.
 * The table also maintains size-binned lists for DP space construction.
 */
public class ClusterTable {

    /**
     * Allocation-free lookup for {@code A \ B}.  A lookup owns one mutable
     * probe key and must therefore stay confined to one worker thread.
     * Successful lookups return the canonical immutable hash already in X.
     */
    public final class ResidualLookup {
        private final ProbeKey probe = new ProbeKey(m);

        private ResidualLookup() {}

        public ClusterHash find(ClusterHash a, ClusterHash b) {
            Entry entry = table.get(probe.setResidual(a, b));
            return entry == null ? null : entry.hash;
        }
    }

    /** Mutable lookup-only key; it is never inserted into {@link #table}. */
    private static final class ProbeKey {
        private final long[] sums;
        private final long[] xors;
        private int size;
        private int hashCode;

        ProbeKey(int m) {
            sums = new long[m];
            xors = new long[m];
        }

        ProbeKey setResidual(ClusterHash a, ClusterHash b) {
            size = a.size - b.size;
            int h = 1;
            for (int s = 0; s < sums.length; s++) {
                long value = a.sums[s] - b.sums[s];
                sums[s] = value;
                h = 31 * h + Long.hashCode(value);
            }
            for (int s = 0; s < xors.length; s++) {
                long value = a.xors[s] ^ b.xors[s];
                xors[s] = value;
                h = 31 * h + Long.hashCode(value);
            }
            hashCode = h;
            return this;
        }

        @Override public int hashCode() { return hashCode; }

        @Override public boolean equals(Object other) {
            if (!(other instanceof ClusterHash cluster) || size != cluster.size) {
                return false;
            }
            return Arrays.equals(sums, cluster.sums)
                && Arrays.equals(xors, cluster.xors);
        }
    }

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

    // The {anchor} singleton hash, for the anchored-outgroup root split
    // null if the anchor
    // taxon appears in no input tree (cannot happen for a real taxon).
    private ClusterHash anchorHash;

    private final int m; // number of hash seeds

    // Anchored-outgroup mode: when
    // true, register only the ANCHOR-FREE orientation of every bipartition (the side
    // not containing the anchor taxon), halving X. Exact when combined with the
    // corresponding anchored DP root.
    private final boolean anchorFreeX;
    private final int     anchor;   // anchor taxon global id (valid iff anchorFreeX)

    // True once any multi-range exemplar has been inserted (consensus emission
    // bridge). Used to gate the GPU weight path, which is single-range-only until
    // the two-tier range-CSR lands.
    private boolean hasMultiRange = false;

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    /** Full registration (both orientations) — backward-compatible entry point. */
    public ClusterTable(List<Tree> trees, PrefixHashArrays pref, int numTaxa) {
        this(trees, pref, numTaxa, false);
    }

    /**
     * @param anchorFreeX  when true, register only the anchor-free orientation of
     *                     every bipartition (halves X). The caller must combine this
     *                     with the anchored DP root.
     */
    public ClusterTable(List<Tree> trees, PrefixHashArrays pref, int numTaxa, boolean anchorFreeX) {
        long t0 = System.nanoTime();
        this.m = pref.numSeeds();
        if (anchorFreeX) {
            throw new IllegalArgumentException(
                "Outgroup anchoring is an unrooted reduction and is not valid in STELAR-X");
        }
        this.anchorFreeX = false;
        this.anchor = -1;

        // Build all-taxa hash from prefix arrays (using any complete tree)
        long[] atSums = new long[m], atXors = new long[m];
        for (int s = 0; s < m; s++) {
            atSums[s] = pref.allTaxaSum(s);
            atXors[s] = pref.allTaxaXor(s);
        }
        allTaxaHash = new ClusterHash(atSums, atXors, numTaxa, m);

        int totalCandidates = 0;
        int treesDone = 0;
        ProgressBar bar = new ProgressBar("Cluster extraction", trees.size());
        for (Tree tree : trees) {
            totalCandidates += extractFromTree(tree, pref, numTaxa);
            bar.update(++treesDone);
        }
        bar.done();

        // Compute the {anchor} singleton hash for the anchored-outgroup root split.
        // Cheap: one range hash from the first tree containing the anchor taxon.
        this.anchorHash = computeAnchorHash(trees, pref, numTaxa);

        // In anchor-free mode the {anchor} singleton (a with-anchor cluster) is NOT
        // registered by walkNodes, but the anchored root's {anchor} child and
        // buildNewick's leaf-name lookup need it — register it explicitly.
        if (anchorFreeX && anchorHash != null) registerAnchorSingleton(trees);

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Cluster extraction: %d candidates -> %d unique clusters in %d ms%s",
            totalCandidates, table.size(), ms, anchorFreeX ? "  [anchor-free X]" : "");
        Logging.debug("  all-taxa cluster (DP root): %s", allTaxaHash);
        if (Logging.isDebug()) {
            logSizeSummary();
        }
    }

    /**
     * Extract clusters from one tree: subtree ranges + complements for every
     * non-root internal node. Leaves are also included (size-1 clusters).
     * Returns the number of candidates generated (before dedup).
     */
    private int extractFromTree(Tree tree, PrefixHashArrays pref, int numTaxa) {
        int ti = tree.treeIndex;
        int L  = tree.leafCount;
        int[] count = {0};
        walkNodes(tree.root, ti, L, -1, pref, numTaxa, count);
        return count[0];
    }

    /**
     * Post-order walk. For every node (including leaves, but excluding root):
     *   full mode        — register the subtree cluster AND its super-complement;
     *   anchor-free mode  — register ONLY the side not containing the anchor taxon.
     */
    private void walkNodes(TreeNode node, int ti, int L, int anchorPos,
                           PrefixHashArrays pref, int numTaxa, int[] count) {
        if (!node.isLeaf()) {
            if (node.isPolytomous()) {
                // Recurse into ALL children. A polytomous node still contributes only
                // its own sub(u) + complement below — NO combo clusters (sub(cᵢ)∪sub(cⱼ))
                // are added; confirmed ASTRAL-MP behaviour (polytomy design §3.3).
                for (TreeNode child : node.children) walkNodes(child, ti, L, anchorPos, pref, numTaxa, count);
            } else {
                walkNodes(node.left,  ti, L, anchorPos, pref, numTaxa, count);
                walkNodes(node.right, ti, L, anchorPos, pref, numTaxa, count);
            }
        }

        if (node.isRoot()) return;  // skip root -- it is the all-taxa cluster

        int lo = node.rangeStart, hi = node.rangeEnd;
        int rangeSize = hi - lo;
        registerCluster(ti, lo, hi, false, rangeSize, L, pref, numTaxa);
        count[0]++;
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
                // Super-complement w.r.t. ALL taxa (S \ [lo,hi))
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

    /**
     * Hash of the {anchor} singleton cluster — the taxon set {@code {anchorTaxon}} —
     * computed as a size-1 range hash from the first input tree that contains the
     * anchor taxon.  Content-based, so it matches the anchor singleton regardless of
     * which tree produced it.  Returns null if the anchor appears in no tree (which
     * cannot happen for a real taxon) or the id is out of range.
     */
    private ClusterHash computeAnchorHash(List<Tree> trees, PrefixHashArrays pref, int numTaxa) {
        int anchor = Config.getInstance().getAnchorTaxon();
        if (anchor < 0 || anchor >= numTaxa) return null;
        for (Tree t : trees) {
            int p = (anchor < t.positionMap.length) ? t.positionMap[anchor] : -1;
            if (p < 0) continue;                       // anchor absent from this tree
            long[] rawSums = new long[m], rawXors = new long[m];
            for (int s = 0; s < m; s++) {
                rawSums[s] = pref.rangeSum(t.treeIndex, s, p, p + 1);
                rawXors[s] = pref.rangeXor(t.treeIndex, s, p, p + 1);
            }
            return new ClusterHash(rawSums, rawXors, 1, m);
        }
        return null;
    }

    /** Register the {anchor} singleton (with-anchor, so not added by walkNodes in anchor-free mode). */
    private void registerAnchorSingleton(List<Tree> trees) {
        for (Tree t : trees) {
            int p = (anchor < t.positionMap.length) ? t.positionMap[anchor] : -1;
            if (p < 0) continue;                       // anchor absent from this tree
            addCluster(anchorHash, new Cluster(t.treeIndex, p, p + 1, false, t.leafCount));
            return;
        }
    }

    // -------------------------------------------------------------------------
    // Queries
    // -------------------------------------------------------------------------

    /** The {anchor} singleton hash for the anchored-outgroup root split (may be null). */
    public ClusterHash getAnchorHash()      { return anchorHash; }
    public Entry get(ClusterHash hash)      { return table.get(hash); }
    public boolean contains(ClusterHash h)  { return table.containsKey(h); }
    public int size()                       { return table.size(); }
    public ClusterHash getAllTaxaHash()     { return allTaxaHash; }
    public Collection<Entry> entries()     { return table.values(); }
    public int numSeeds()                  { return m; }

    /** Create one allocation-free residual lookup for a single worker thread. */
    public ResidualLookup newResidualLookup() { return new ResidualLookup(); }

    /** All cluster hashes of a given size. */
    public List<ClusterHash> getBySize(int sz) {
        return sizeBins.getOrDefault(sz, Collections.emptyList());
    }

    /** All sizes present in X. */
    public Set<Integer> sizes() { return sizeBins.keySet(); }

    /**
     * Add all bipartitions from one additional tree (e.g. the UPGMA guide tree)
     * into the cluster table.  The root's bipartition (all-taxa cluster) is
     * automatically skipped, exactly as in the constructor's per-tree walk.
     */
    public void addTree(Tree tree, PrefixHashArrays pref, int numTaxa) {
        extractFromTree(tree, pref, numTaxa);
    }

    /**
     * Insert a cluster whose {@link ClusterHash} is already known — used by the
     * consensus polytomy-resolution emission bridge, where the signature was
     * computed via the consensus tree's own prefix scans (and is comparable to
     * gene-tree signatures because the same {@code TaxonHasher} produced both).
     *
     * Mirrors {@link #registerCluster}'s table + sizeBins bookkeeping so the new
     * cluster participates in Mode 2 cross-tree transitions exactly like a
     * gene-tree-derived cluster. The all-taxa cluster and existing hashes are
     * skipped. The exemplar may be single- or multi-range.
     *
     * @return true iff newly added
     */
    public boolean addCluster(ClusterHash hash, Cluster exemplar) {
        if (hash.equals(allTaxaHash)) return false;
        if (table.containsKey(hash))  return false;
        table.put(hash, new Entry(hash, exemplar));
        sizeBins.computeIfAbsent(hash.size, k -> new ArrayList<>()).add(hash);
        if (exemplar.isMultiRange()) hasMultiRange = true;
        return true;
    }

    /** True iff any inserted exemplar is multi-range (gates the GPU weight path). */
    public boolean hasMultiRange() { return hasMultiRange; }

    // -------------------------------------------------------------------------

    private void logSizeSummary() {
        if (sizeBins.isEmpty()) return;
        int minSz = Integer.MAX_VALUE, maxSz = 0;
        for (int sz : sizeBins.keySet()) {
            minSz = Math.min(minSz, sz);
            maxSz = Math.max(maxSz, sz);
        }
        Logging.debug("  cluster size range: [%d, %d]", minSz, maxSz);
        Logging.debug("  distinct sizes: %d", sizeBins.size());
    }
}
