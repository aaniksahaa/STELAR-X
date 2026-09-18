package stelarx.greedy;

import stelarx.cluster.Cluster;
import stelarx.tree.Tree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * INSERT procedure for the laminar forest — directly mirrors the structure of
 * ASTRAL-MP's {@code Utils.buildTreeFromClusters} loop body:
 *
 *   1. Skip clusters of size ≤ 1 or size == n (trivial: already in the star).
 *   2. Find LCA(C) in the current laminar forest.
 *   3. For every child of LCA, count how many of C's taxa fall in its subtree.
 *   4. A touched child whose hit count != size is a CROSS-CUT → REJECT.
 *   5. If every touched child is wholly covered AND
 *         #touched-children ≥ 2 AND
 *         touched != all of LCA's children
 *      → create a new internal node, adopt the touched children, ACCEPT.
 *      Otherwise (single or all-children touched) → SKIP (no new bipartition).
 *
 * Outcomes are exposed via {@link Outcome} so the oracle and the driver can
 * cross-check or observe accept/reject/skip decisions.
 */
final class LaminarBuilder {

    enum Outcome { ACCEPT, SKIP_TRIVIAL, SKIP_REDUNDANT, REJECT_CROSS_CUT }

    private final LaminarForest forest;
    private final List<Tree> trees;
    private final int numTaxa;

    LaminarBuilder(LaminarForest forest, List<Tree> trees, int numTaxa) {
        this.forest  = forest;
        this.trees   = trees;
        this.numTaxa = numTaxa;
    }

    /**
     * Attempt to refine the laminar forest by inserting bipartition b's
     * canonical-side cluster.
     */
    Outcome insert(Bipartition b) {
        Cluster c = b.canonicalExemplar;

        // (1) Trivial-size shortcut — matches ASTRAL-MP's
        //     `tc.getClusterSize() <= 1 || tc.getClusterSize() == taxonCount()` filter.
        if (b.size <= 1 || b.size == numTaxa) {
            return Outcome.SKIP_TRIVIAL;
        }

        // (2) Enumerate the cluster's taxa from its exemplar tree's postorder.
        int[] taxa = enumerateTaxa(c);

        // (3) LCA of the cluster's leaves in the laminar forest.
        LaminarNode lca = forest.lcaOfTaxa(taxa);

        // (4) Per-child-of-LCA hit counts.
        //     hits maps child node id → count of C's taxa under that child.
        HashMap<LaminarNode, int[]> hits = new HashMap<>(8);
        for (int t : taxa) {
            LaminarNode child = forest.childOnPathTo(lca, t);
            int[] count = hits.get(child);
            if (count == null) {
                count = new int[]{0};
                hits.put(child, count);
            }
            count[0]++;
        }

        // (5) Whole-node test.  Mirrors ASTRAL-MP's
        //     `if (tc.containsCluster(childCluster)) moved += child`
        //     with the post-check `remainingleaves != 0` → reject.
        List<LaminarNode> moved = new ArrayList<>(hits.size());
        int accountedFor = 0;
        for (var entry : hits.entrySet()) {
            LaminarNode child = entry.getKey();
            int hit = entry.getValue()[0];
            if (hit != child.size) {
                // Some leaf in this child is in C, some isn't → cross-cut.
                return Outcome.REJECT_CROSS_CUT;
            }
            moved.add(child);
            accountedFor += hit;
        }
        if (accountedFor != taxa.length) {
            // Defensive: should be impossible if LCA is correct, but assert.
            return Outcome.REJECT_CROSS_CUT;
        }

        // (6) Result classification.
        //
        //   moved.size() == 0       → impossible (every taxon has a child) so
        //                             not handled separately
        //   moved.size() == 1       → C exactly matches an existing child; the
        //                             new node would have one child, no new
        //                             bipartition.  SKIP.
        //   moved == LCA.children   → C == taxa(LCA), redundant w/ LCA itself.
        //   otherwise               → ACCEPT.
        if (moved.size() < 2) {
            return Outcome.SKIP_REDUNDANT;
        }
        if (moved.size() == lca.children.size()) {
            return Outcome.SKIP_REDUNDANT;
        }

        // (7) Create new internal node, adopt moved children.
        LaminarNode w = forest.createInternal(lca);
        lca.children.add(w);
        forest.moveChildren(lca, w, moved);
        return Outcome.ACCEPT;
    }

    /**
     * Enumerate the taxon ids in a cluster's exemplar.
     *
     * For {@code complement=true} the cluster represents the SUPER-complement
     * w.r.t. all n taxa (i.e. {0..n-1} \ {taxa in the in-tree range}).  When
     * the source tree is incomplete this side includes the taxa missing from
     * that tree — we walk taxon ids globally and use {@link Tree#positionMap}
     * to detect "in the range" rather than the in-tree postorder.
     */
    private int[] enumerateTaxa(Cluster ex) {
        Tree tree = trees.get(ex.treeIndex);
        int[] arr = tree.postorderArray;
        int[] out = new int[ex.size];
        int k = 0;
        if (!ex.complement) {
            for (int i = ex.left; i < ex.right; i++) out[k++] = arr[i];
        } else {
            int[] pos = tree.positionMap;
            for (int t = 0; t < numTaxa; t++) {
                int p = pos[t];
                if (p < 0 || p < ex.left || p >= ex.right) out[k++] = t;
            }
        }
        if (k != ex.size) {
            throw new IllegalStateException(
                "enumerateTaxa size mismatch: produced " + k + ", expected " + ex.size);
        }
        return out;
    }
}
