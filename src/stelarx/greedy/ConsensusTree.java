package stelarx.greedy;

import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.NewickLabel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.TreeSet;
import java.util.function.Consumer;

/**
 * Frozen snapshot of a {@link LaminarForest}'s state at a given threshold.
 *
 * Part I provides:
 *   - the {@link SNode} tree itself
 *   - Newick + canonical-leaf-set dumps (for diffing snapshots between the
 *     fast laminar builder and the brute-force oracle)
 *
 * Part II additions (this file):
 *   - the postorder leaf array {@code aCons[0..n-1]} (§7.1)
 *   - prefix-scan hash arrays {@code p1[s][0..n], p2[s][0..n]} per seed s,
 *     applying the SAME taxon hasher used to build the gene-tree prefix arrays
 *     so a consensus-derived signature matches a gene-tree-derived signature
 *     for the same taxon set
 *   - per-node range stamps {@code (lo, hi)} into {@code aCons}
 *   - O(1) range signature queries {@link #sigma1}/{@link #sigma2}, O(1)
 *     complement via the whole-tree totals, and O(m) disjoint multi-range
 *     unions through the (additive, XOR) group operators
 *
 * The cross-source signature parity is the entire reason the polytomy
 * resolution can dedupe its emitted bipartitions against X without any
 * encoding conversion.
 */
public final class ConsensusTree {

    /**
     * One node of the snapshot.  Each node has a dense integer id, a taxonId
     * (≥ 0 for leaves, -1 for internal nodes), and a {@code [lo, hi)} range
     * into {@code aCons} stamped at construction time.
     */
    public static final class SNode {
        public final int id;
        public final int taxonId;            // -1 if internal
        public final List<SNode> children;   // empty for leaves
        /** Range start (inclusive) into the parent ConsensusTree's aCons. */
        int lo;
        /** Range end (exclusive) into the parent ConsensusTree's aCons. */
        int hi;

        SNode(int id, int taxonId, List<SNode> children) {
            this.id = id;
            this.taxonId = taxonId;
            this.children = children;
        }

        public boolean isLeaf()    { return taxonId >= 0; }
        public int rangeLo()       { return lo; }
        public int rangeHi()       { return hi; }
        public int rangeSize()     { return hi - lo; }
    }

    // ── Tree structure ──
    private final SNode root;
    private final int numTaxa;
    private final int numNodes;
    private final int numInternalNodes;     // excluding virtual root
    private final int numPolytomies;        // internal nodes with > 2 children

    // ── Prefix-scan hashing over consensus postorder ──
    private final int[] aCons;              // [n] leaf id at each postorder position
    private final int   m;                  // number of hash seeds
    private final long[][] p1;              // [m][n+1] prefix sum, p1[s][0] = 0
    private final long[][] p2;              // [m][n+1] prefix XOR, p2[s][0] = 0

    private ConsensusTree(SNode root, int numTaxa, int numNodes,
                          int numInternalNodes, int numPolytomies,
                          int[] aCons, long[][] p1, long[][] p2) {
        this.root             = root;
        this.numTaxa          = numTaxa;
        this.numNodes         = numNodes;
        this.numInternalNodes = numInternalNodes;
        this.numPolytomies    = numPolytomies;
        this.aCons            = aCons;
        this.p1               = p1;
        this.p2               = p2;
        this.m                = p1.length;
    }

    // ── Construction ───────────────────────────────────────────────────────

    /**
     * Take a snapshot of the current laminar-forest state.  Deep-copies the
     * children structure so subsequent INSERTs into the source forest do not
     * mutate the snapshot, then performs a single postorder pass to compute
     * {@code aCons[]} and stamp every SNode's {@code (lo, hi)} range.
     *
     * @param hasher must be the same TaxonHasher that produced the gene-tree
     *               prefix arrays — otherwise cross-source signature parity
     *               is broken.
     */
    public static ConsensusTree snapshot(LaminarForest forest, TaxonHasher hasher) {
        int n = forest.numTaxa;
        int[] counters = new int[3];   // [0]=nextId, [1]=internals, [2]=polytomies
        IdentityHashMap<LaminarNode, Integer> minTaxon = new IdentityHashMap<>();
        computeMinTaxon(forest.virtualRoot, minTaxon);
        SNode root = copyRec(forest.virtualRoot, counters, minTaxon);

        int[] aCons = new int[n];
        int[] pos = {0};
        postOrderStamp(root, aCons, pos);
        if (pos[0] != n) {
            throw new IllegalStateException(
                "ConsensusTree.snapshot: aCons fill " + pos[0] + " != n=" + n);
        }

        int m = hasher.numSeeds();
        long[][] p1 = new long[m][n + 1];
        long[][] p2 = new long[m][n + 1];
        for (int s = 0; s < m; s++) {
            long sum = 0L, xor = 0L;
            for (int i = 0; i < n; i++) {
                long h = hasher.get(s, aCons[i]);
                sum += h;       // unsigned add mod 2^64 (Java long wraps)
                xor ^= h;
                p1[s][i + 1] = sum;
                p2[s][i + 1] = xor;
            }
        }

        return new ConsensusTree(root, n, counters[0], counters[1], counters[2],
                                  aCons, p1, p2);
    }

    /** Precompute a canonical, order-independent key for every disjoint child. */
    private static int computeMinTaxon(LaminarNode src,
                                       IdentityHashMap<LaminarNode, Integer> minima) {
        int min = src.isLeaf() ? src.taxonId : Integer.MAX_VALUE;
        for (LaminarNode child : src.children) {
            min = Math.min(min, computeMinTaxon(child, minima));
        }
        minima.put(src, min);
        return min;
    }

    private static SNode copyRec(LaminarNode src, int[] counters,
                                 IdentityHashMap<LaminarNode, Integer> minima) {
        int id = counters[0]++;
        if (src.isLeaf()) {
            return new SNode(id, src.taxonId, Collections.emptyList());
        }
        // Laminar insertion may move children in thread-dependent encounter order.
        // The children are disjoint, so their minimum taxon IDs are unique and give
        // a cheap canonical ordering.  This stabilises node IDs, aCons ranges, and
        // the mapping from deterministic RNG draws to polytomy groups.
        List<LaminarNode> ordered = new ArrayList<>(src.children);
        ordered.sort(Comparator.comparingInt(minima::get));
        List<SNode> kids = new ArrayList<>(ordered.size());
        for (LaminarNode c : ordered) kids.add(copyRec(c, counters, minima));
        // Don't count the virtual root itself
        if (src.parent >= 0) {
            counters[1]++;
            if (kids.size() > 2) counters[2]++;
        }
        return new SNode(id, -1, kids);
    }

    private static void postOrderStamp(SNode n, int[] aCons, int[] pos) {
        int lo = pos[0];
        if (n.isLeaf()) {
            aCons[pos[0]++] = n.taxonId;
        } else {
            for (SNode c : n.children) postOrderStamp(c, aCons, pos);
        }
        n.lo = lo;
        n.hi = pos[0];
    }

    // ── Accessors ──────────────────────────────────────────────────────────

    public int numTaxa()           { return numTaxa; }
    public int numNodes()          { return numNodes; }
    public int numInternalNodes()  { return numInternalNodes; }
    public int numPolytomies()     { return numPolytomies; }
    public int numSeeds()          { return m; }
    public SNode root()            { return root; }
    public int[] aCons()           { return aCons; }

    /** Visit every non-leaf node (including the virtual root) in preorder. */
    public void forEachInternalNode(Consumer<SNode> consumer) {
        visitInternal(root, consumer);
    }

    private void visitInternal(SNode n, Consumer<SNode> c) {
        if (n.isLeaf()) return;
        c.accept(n);
        for (SNode child : n.children) visitInternal(child, c);
    }

    // ── O(1) signature queries ─────────────────────────────────────────────

    /** ϕ1 signature of the contiguous range {@code aCons[lo..hi)}. */
    public long sigma1(int s, int lo, int hi) {
        return p1[s][hi] - p1[s][lo];
    }

    /** ϕ2 signature of the contiguous range {@code aCons[lo..hi)}. */
    public long sigma2(int s, int lo, int hi) {
        return p2[s][hi] ^ p2[s][lo];
    }

    /** ϕ1 signature of the WHOLE tree (= sum of all taxon hashes). */
    public long totalSigma1(int s) { return p1[s][numTaxa]; }
    /** ϕ2 signature of the WHOLE tree. */
    public long totalSigma2(int s) { return p2[s][numTaxa]; }

    /** ϕ1 signature of the complement of {@code aCons[lo..hi)}. */
    public long complementSigma1(int s, int lo, int hi) {
        return p1[s][numTaxa] - sigma1(s, lo, hi);
    }
    /** ϕ2 signature of the complement of {@code aCons[lo..hi)}. */
    public long complementSigma2(int s, int lo, int hi) {
        return p2[s][numTaxa] ^ sigma2(s, lo, hi);
    }

    /**
     * Combine pairwise-disjoint ranges into one signature.  Required for
     * resolutions whose side spans several non-adjacent groups (Step B's
     * induced split, Step A's union of non-adjacent groups).
     *
     * PRECONDITION: the ranges must be pairwise disjoint.  ϕ1 additive
     * double-counts overlaps; ϕ2 XOR silently cancels duplicated taxa.
     * The caller MUST guarantee this — which is automatic for unions of
     * children of a single polytomy, since children's ranges tile [v.lo, v.hi)
     * without overlap, but is NOT automatic for arbitrary callers.
     */
    public long combineDisjointSigma1(int s, int[] los, int[] his) {
        long sum = 0L;
        for (int j = 0; j < los.length; j++) {
            sum += p1[s][his[j]] - p1[s][los[j]];
        }
        return sum;
    }

    /** XOR counterpart of {@link #combineDisjointSigma1}. */
    public long combineDisjointSigma2(int s, int[] los, int[] his) {
        long xor = 0L;
        for (int j = 0; j < los.length; j++) {
            xor ^= p2[s][his[j]] ^ p2[s][los[j]];
        }
        return xor;
    }

    // ── Existing dump APIs (Newick + canonical leaf sets) ──────────────────

    /** Newick string using taxon names; no branch lengths. */
    public String toNewick(TaxonRegistry registry) {
        StringBuilder sb = new StringBuilder();
        writeNewick(root, registry, sb);
        sb.append(';');
        return sb.toString();
    }

    private void writeNewick(SNode n, TaxonRegistry reg, StringBuilder sb) {
        if (n.isLeaf()) {
            sb.append(NewickLabel.encode(reg.getName(n.taxonId)));
            return;
        }
        sb.append('(');
        for (int i = 0; i < n.children.size(); i++) {
            if (i > 0) sb.append(',');
            writeNewick(n.children.get(i), reg, sb);
        }
        sb.append(')');
    }

    /**
     * "{1,4,7}\n{2,5}\n..." dump — one line per non-root internal node, leaf
     * ids sorted, lines sorted lexicographically.
     */
    public String canonicalLeafSets() {
        List<String> lines = new ArrayList<>();
        TreeSet<Integer> scratch = new TreeSet<>();
        collectLeafSets(root, /*isRoot=*/true, lines, scratch);
        Collections.sort(lines);
        return String.join("\n", lines);
    }

    private void collectLeafSets(SNode n, boolean isRoot,
                                  List<String> out, TreeSet<Integer> scratch) {
        if (n.isLeaf()) return;
        if (!isRoot) {
            scratch.clear();
            for (int p = n.lo; p < n.hi; p++) scratch.add(aCons[p]);
            StringBuilder sb = new StringBuilder("{");
            boolean first = true;
            for (int t : scratch) {
                if (!first) sb.append(',');
                sb.append(t);
                first = false;
            }
            sb.append('}');
            out.add(sb.toString());
        }
        for (SNode c : n.children) collectLeafSets(c, false, out, scratch);
    }
}
