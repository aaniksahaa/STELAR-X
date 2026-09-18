package stelarx.partition;

import stelarx.Logging;
import stelarx.cluster.ClusterHash;
import stelarx.util.ProgressBar;
import stelarx.hash.PrefixHashArrays;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.*;

/**
 * Table of unique rooted gene-tree child partitions with their frequencies.
 *
 * For each internal node u of each rooted gene tree g (root included) we extract:
 *   children(u) = the unordered collection of its child-subtree leaf sets.
 *
 * Deduplication: PartitionHash is order-invariant over the actual children and
 * excludes all ambient outside taxa. Nodes that cannot display a resolved rooted
 * triplet are omitted entirely.
 */
public class PartitionTable {

    public static final class Entry {
        public final PartitionHash hash;
        public final Partition     exemplar;
        public int                 frequency;

        Entry(PartitionHash h, Partition p) { this.hash = h; this.exemplar = p; this.frequency = 1; }
    }

    private final Map<PartitionHash, Entry> table = new HashMap<>();
    private final int m;
    private boolean hasPoly = false;   // true once any d>3 (polytomous) partition is stored

    /** True iff any extracted partition is polytomous (d > 3). */
    public boolean hasPolytomousPartitions() { return hasPoly; }

    // -------------------------------------------------------------------------

    public PartitionTable(List<Tree> trees, PrefixHashArrays pref) {
        long t0 = System.nanoTime();
        this.m = pref.numSeeds();

        int totalCandidates = 0;
        int zeroCandidates = 0;
        int treesDone = 0;
        ProgressBar bar = new ProgressBar("Child-partition extraction", trees.size());
        for (Tree tree : trees) {
            int[] counts = extractFromTree(tree, pref);
            totalCandidates += counts[0];
            zeroCandidates += counts[1];
            bar.update(++treesDone);
        }
        bar.done();

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Rooted partition extraction: %d internal nodes -> %d contributing occurrences "
                + "(%d identically zero skipped) -> %d unique child partitions in %d ms",
            totalCandidates, totalCandidates - zeroCandidates, zeroCandidates, table.size(), ms);
    }

    /** Returns {all internal nodes, identically-zero nodes skipped}. */
    private int[] extractFromTree(Tree tree, PrefixHashArrays pref) {
        int ti = tree.treeIndex;
        int[] counts = {0, 0};
        extractNode(tree.root, ti, pref, counts);
        return counts;
    }

    /**
     * Recurse post-order. For every internal node u, root included, register the
     * unordered collection of actual child-subtree sets.
     */
    private void extractNode(TreeNode node, int ti, PrefixHashArrays pref, int[] counts) {
        if (node.isLeaf()) return;
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) extractNode(child, ti, pref, counts);
        } else {
            extractNode(node.left,  ti, pref, counts);
            extractNode(node.right, ti, pref, counts);
        }
        counts[0]++;

        // ── Polytomous node: unordered collection of k child subtrees ──
        if (node.isPolytomous()) {
            int k = node.children.length;
            ClusterHash[] hashes = new ClusterHash[k];
            int[] sizes      = new int[k];
            int[] partStarts = new int[k];
            int[] partEnds   = new int[k];
            for (int i = 0; i < k; i++) {
                TreeNode c = node.children[i];
                int cs = c.rangeStart, ce = c.rangeEnd, szi = ce - cs;
                partStarts[i] = cs; partEnds[i] = ce; sizes[i] = szi;
                hashes[i] = buildHash(ti, cs, ce, szi, pref);
            }
            boolean canDisplayTriplet = false;
            for (int i = 0; i < k; i++) {
                if (sizes[i] >= 2) { canDisplayTriplet = true; break; }
            }
            if (!canDisplayTriplet) {
                counts[1]++;
                return;
            }
            PartitionHash ph = new PartitionHash(hashes);
            Entry existing = table.get(ph);
            if (existing != null) {
                existing.frequency++;
            } else {
                Partition p = new Partition(sizes, partStarts, partEnds, ti);
                table.put(ph, new Entry(ph, p));
                hasPoly = true;
            }
            return;
        }

        // ── Binary node: unordered pair of actual child subtrees ──
        int lStart = node.left.rangeStart,  lEnd = node.left.rangeEnd;
        int rStart = node.right.rangeStart, rEnd = node.right.rangeEnd;
        int sz1 = lEnd - lStart;
        int sz2 = rEnd - rStart;

        // Two singleton children cannot supply the paired leaves of any resolved
        // rooted triplet, so this node contributes zero for every candidate tree.
        if (sz1 == 1 && sz2 == 1) {
            counts[1]++;
            return;
        }

        ClusterHash h1 = buildHash(ti, lStart, lEnd, sz1, pref);
        ClusterHash h2 = buildHash(ti, rStart, rEnd, sz2, pref);

        PartitionHash ph = new PartitionHash(h1, h2);

        Entry existing = table.get(ph);
        if (existing != null) {
            existing.frequency++;
        } else {
            Partition p = new Partition(sz1, sz2, ti, lStart, lEnd, rStart, rEnd);
            table.put(ph, new Entry(ph, p));
        }
    }

    /** Compute a ClusterHash for one child-subtree range. */
    private ClusterHash buildHash(int ti, int lo, int hi, int size, PrefixHashArrays pref) {
        long[] rawSums = new long[m], rawXors = new long[m];
        for (int s = 0; s < m; s++) {
            rawSums[s] = pref.rangeSum(ti, s, lo, hi);
            rawXors[s] = pref.rangeXor(ti, s, lo, hi);
        }
        return new ClusterHash(rawSums, rawXors, size, m);
    }

    // -------------------------------------------------------------------------

    public Entry get(PartitionHash ph) { return table.get(ph); }
    public int size()                  { return table.size(); }
    public Collection<Entry> entries() { return table.values(); }
}
