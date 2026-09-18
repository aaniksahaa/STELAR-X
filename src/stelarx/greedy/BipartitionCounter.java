package stelarx.greedy;

import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

/**
 * Phase 1 of greedy consensus: collect unique unrooted bipartitions and sort
 * them by frequency descending.
 *
 * Two collection sources are supported:
 *
 *   {@link #collectFromGeneTrees}   — walks the supplied gene trees directly,
 *                                     emitting one cluster per internal node and
 *                                     deduping (cluster vs. its complement) into
 *                                     a fresh per-bipartition count.  This is
 *                                     the head-to-head-faithful path: it matches
 *                                     ASTRAL-MP's {@code Utils.getGeneClusters}
 *                                     + count HashMap, and excludes anything
 *                                     not in the gene trees (e.g. the UPGMA
 *                                     guide tree).
 *
 *   {@link #collectFromClusterTable} — aggregates the two-sides-per-bipartition
 *                                      frequencies that ClusterTable already
 *                                      maintains.  This is the fast path when
 *                                      the caller is OK with whatever ClusterTable
 *                                      was built from (gene trees + UPGMA).
 *
 * Both return a sorted list of {@link Bipartition}s with the canonical-side
 * selection rule:
 *   - If both sides have size in [2, n-1]: pick smaller; tie-break by hash.
 *   - If one side has size 1: pick the OTHER side (size n-1).  ASTRAL-MP's
 *     buildTreeFromClusters skips clusters of size ≤ 1, so picking the size-1
 *     side would silently drop this bipartition.
 *   - If both sides have size 1 (n == 2): no bipartition can be formed; skip.
 *
 * The list is sorted by frequency descending; ties are broken by canonical
 * hash for determinism.
 */
public final class BipartitionCounter {

    private BipartitionCounter() {}

    // ── Mutable accumulator used by collectFromGeneTrees ─────────────────────
    private static final class Accum {
        final ClusterHash subHash;   // the (subtree-side) hash registered first
        final Cluster     exemplar;  // pointing at the subtree-side range
        final int         subSize;   // |subtree side|
        int               freq;

        Accum(ClusterHash subHash, Cluster exemplar, int subSize) {
            this.subHash = subHash;
            this.exemplar = exemplar;
            this.subSize = subSize;
            this.freq = 1;
        }
    }

    /**
     * Walk the supplied gene trees directly and produce per-bipartition counts.
     * This is the ASTRAL-MP-faithful path — used to exclude the UPGMA guide
     * tree's contributions from the greedy-consensus frequencies.
     *
     * @param geneTrees    gene trees only (no UPGMA, no padding)
     * @param pref         prefix hash arrays — must cover the supplied trees'
     *                     {@link Tree#treeIndex} values
     * @param allTaxaHash  the all-taxa cluster signature (typically
     *                     {@code ClusterTable.getAllTaxaHash()})
     * @param numTaxa      n
     */
    public static List<Bipartition> collectFromGeneTrees(
            List<Tree> geneTrees, PrefixHashArrays pref,
            ClusterHash allTaxaHash, int numTaxa) {

        if (numTaxa < 2) return new ArrayList<>(0);
        int m = pref.numSeeds();

        HashMap<ClusterHash, Accum> count = new HashMap<>(geneTrees.size() * numTaxa);

        for (Tree tree : geneTrees) {
            walkInternals(tree, count, pref, allTaxaHash, m, numTaxa);
        }

        return toSortedList(count, allTaxaHash, numTaxa);
    }

    private static void walkInternals(Tree tree, HashMap<ClusterHash, Accum> count,
                                       PrefixHashArrays pref, ClusterHash allTaxaHash,
                                       int m, int numTaxa) {
        walkRec(tree.root, tree, count, pref, allTaxaHash, m, numTaxa);
    }

    private static void walkRec(TreeNode node, Tree tree,
                                 HashMap<ClusterHash, Accum> count,
                                 PrefixHashArrays pref, ClusterHash allTaxaHash,
                                 int m, int numTaxa) {
        if (node.isLeaf()) return;
        walkRec(node.left,  tree, count, pref, allTaxaHash, m, numTaxa);
        walkRec(node.right, tree, count, pref, allTaxaHash, m, numTaxa);
        if (node.isRoot()) return;

        int lo = node.rangeStart, hi = node.rangeEnd;
        int subSize = hi - lo;
        // Internal-node subtree size is always in [2, leafCount−1].  For the
        // ASTRAL-MP-faithful semantics we treat the cluster as the in-tree
        // range; bipartitions where the complement lies entirely in missing
        // taxa (subSize == leafCount) are still legitimate w.r.t. all n taxa.
        if (subSize < 2) return;
        if (subSize >= numTaxa) return;

        // Compute the subtree-side hash via prefix arrays.
        long[] rawSums = new long[m], rawXors = new long[m];
        int ti = tree.treeIndex;
        for (int s = 0; s < m; s++) {
            rawSums[s] = pref.rangeSum(ti, s, lo, hi);
            rawXors[s] = pref.rangeXor(ti, s, lo, hi);
        }
        ClusterHash subHash = new ClusterHash(rawSums, rawXors, subSize, m);

        // Dedupe against existing entries and their complements.
        Accum existing = count.get(subHash);
        if (existing != null) {
            existing.freq++;
            return;
        }
        ClusterHash compHash = ClusterHash.residual(allTaxaHash, subHash);
        Accum compExisting = count.get(compHash);
        if (compExisting != null) {
            compExisting.freq++;
            return;
        }

        // New bipartition.  Exemplar points at the in-tree range (complement=false).
        Cluster ex = new Cluster(ti, lo, hi, false, tree.leafCount);
        count.put(subHash, new Accum(subHash, ex, subSize));
    }

    /**
     * Aggregate existing per-side frequencies from {@link ClusterTable} into
     * per-bipartition entries.  Includes whatever ClusterTable was built from
     * (gene trees + UPGMA guide tree, in the standard STELAR-X pipeline).
     */
    public static List<Bipartition> collectFromClusterTable(
            ClusterTable clusterTable, int numTaxa) {

        if (numTaxa < 2) return new ArrayList<>(0);
        ClusterHash allTaxa = clusterTable.getAllTaxaHash();

        List<Bipartition> out = new ArrayList<>(clusterTable.size() / 2 + 16);
        HashSet<ClusterHash> consumed = new HashSet<>(clusterTable.size() * 2);

        for (ClusterTable.Entry e : clusterTable.entries()) {
            if (consumed.contains(e.hash)) continue;
            int sizeA = e.hash.size;
            int sizeB = numTaxa - sizeA;

            ClusterHash compHash = ClusterHash.residual(allTaxa, e.hash);
            ClusterTable.Entry compEntry = clusterTable.get(compHash);
            consumed.add(e.hash);
            if (compEntry != null) consumed.add(compEntry.hash);

            ClusterTable.Entry chosen;
            int chosenSize;
            if (sizeA == 1 && sizeB == 1) continue;
            else if (sizeA == 1) {
                if (compEntry == null) continue;
                chosen = compEntry; chosenSize = sizeB;
            } else if (sizeB == 1) {
                chosen = e; chosenSize = sizeA;
            } else {
                if (compEntry == null || sizeA < sizeB) {
                    chosen = e; chosenSize = sizeA;
                } else if (sizeA > sizeB) {
                    chosen = compEntry; chosenSize = sizeB;
                } else {
                    chosen = (compareHash(e.hash, compEntry.hash) <= 0) ? e : compEntry;
                    chosenSize = sizeA;
                }
            }
            int freq = (compEntry != null)
                ? Math.max(e.frequency, compEntry.frequency)
                : e.frequency;
            out.add(new Bipartition(chosen.hash, chosen.exemplar, chosenSize, freq));
        }

        out.sort((a, b) -> {
            int byFreq = Integer.compare(b.frequency, a.frequency);
            if (byFreq != 0) return byFreq;
            return compareHash(a.canonicalHash, b.canonicalHash);
        });
        return out;
    }

    /** Convert the gene-tree-walk accumulator into a sorted Bipartition list. */
    private static List<Bipartition> toSortedList(
            HashMap<ClusterHash, Accum> count, ClusterHash allTaxaHash, int numTaxa) {

        List<Bipartition> out = new ArrayList<>(count.size());
        for (Accum a : count.values()) {
            int sizeA = a.subSize;
            int sizeB = numTaxa - sizeA;
            if (sizeA == 1 && sizeB == 1) continue;

            ClusterHash canonHash;
            Cluster canonExemplar;
            int canonSize;

            ClusterHash compHash = ClusterHash.residual(allTaxaHash, a.subHash);

            // Determine canonical side
            boolean useA;
            if (sizeA == 1)         useA = false;   // complement is size n-1
            else if (sizeB == 1)    useA = true;
            else if (sizeA < sizeB) useA = true;
            else if (sizeA > sizeB) useA = false;
            else                    useA = (compareHash(a.subHash, compHash) <= 0);

            if (useA) {
                canonHash = a.subHash;
                canonExemplar = a.exemplar;
                canonSize = sizeA;
            } else {
                canonHash = compHash;
                canonExemplar = new Cluster(a.exemplar.treeIndex,
                                            a.exemplar.left, a.exemplar.right,
                                            /*complement=*/true, numTaxa);
                canonSize = sizeB;
            }
            out.add(new Bipartition(canonHash, canonExemplar, canonSize, a.freq));
        }

        out.sort((x, y) -> {
            int byFreq = Integer.compare(y.frequency, x.frequency);
            if (byFreq != 0) return byFreq;
            return compareHash(x.canonicalHash, y.canonicalHash);
        });
        return out;
    }

    /** Total deterministic compare over (sums[], xors[], size). */
    private static int compareHash(ClusterHash a, ClusterHash b) {
        int n1 = a.sums.length, n2 = b.sums.length;
        int n = Math.min(n1, n2);
        for (int s = 0; s < n; s++) {
            int c = Long.compareUnsigned(a.sums[s], b.sums[s]);
            if (c != 0) return c;
        }
        for (int s = 0; s < n; s++) {
            int c = Long.compareUnsigned(a.xors[s], b.xors[s]);
            if (c != 0) return c;
        }
        int c = Integer.compare(n1, n2);
        if (c != 0) return c;
        return Integer.compare(a.size, b.size);
    }
}
