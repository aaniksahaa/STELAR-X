package stelarx.weight;

import stelarx.tree.Tree;

/**
 * O(min(|P|, |Q|)) intersection counts between postorder ranges.
 *
 * For complete gene trees (Lg = L), complement intersections reduce to:
 *   |comp(P) ∩ Q| = |Q| - |P ∩ Q|
 *   |P ∩ comp(Q)| = |P| - |P ∩ Q|
 *   |comp(P) ∩ comp(Q)| = n - |P| - |Q| + |P ∩ Q|
 *
 * So we only ever need the one "core" non-complement count.
 */
public final class IntersectionCounter {

    private IntersectionCounter() {}

    /**
     * |range_in_treeA ∩ range_in_treeB| -- both ranges non-complement.
     * Iterates over the smaller range; looks up each taxon in the other tree's
     * positionMap to check membership.
     */
    public static int coreIntersect(Tree tA, int loA, int hiA,
                                     Tree tB, int loB, int hiB) {
        int count = 0;
        if (hiA - loA <= hiB - loB) {
            for (int pos = loA; pos < hiA; pos++) {
                int taxon = tA.postorderArray[pos];
                int posB  = tB.positionMap[taxon];
                if (posB >= loB && posB < hiB) count++;
            }
        } else {
            for (int pos = loB; pos < hiB; pos++) {
                int taxon = tB.postorderArray[pos];
                int posA  = tA.positionMap[taxon];
                if (posA >= loA && posA < hiA) count++;
            }
        }
        return count;
    }

    /**
     * |M_range ∩ cluster| where M_range is non-complement (a gene-tree subtree
     * range) and cluster may be complement (a species-tree candidate part).
     *
     * For complete trees: |comp(cluster_range) ∩ M_range| = |M_range| - |cluster_range ∩ M_range|
     *
     * @param tGT       gene-tree tree
     * @param loGT, hiGT  range in gene-tree postorder array (non-complement)
     * @param tC        exemplar tree of the candidate cluster
     * @param loC, hiC  range in candidate cluster's exemplar tree
     * @param cComp     whether the candidate cluster is complement w.r.t. its tree
     * @param sizeGTRange  = hiGT - loGT (passed in to avoid recomputation)
     */
    public static int intersect(Tree tGT, int loGT, int hiGT,
                                 Tree tC, int loC, int hiC, boolean cComp,
                                 int sizeGTRange) {
        int core = coreIntersect(tGT, loGT, hiGT, tC, loC, hiC);
        return cComp ? (sizeGTRange - core) : core;
    }

    /**
     * |A ∩ Lg_GT| — how many of the gene-tree's leaves fall inside cluster A.
     *
     * For a super-complement cluster A = S\sub(u):
     *   |A ∩ Lg_GT| = L_GT - |sub(u) ∩ Lg_GT|  (since M ⊆ Lg_GT ⊆ S)
     * For a plain cluster A = sub(u):
     *   |A ∩ Lg_GT| = |sub(u) ∩ Lg_GT|
     *
     * Used to compute the correct row sum for the A-row of the intersection matrix
     * when tGT is an incomplete gene tree (L_GT < n).
     * For complete gene trees (L_GT == n) this returns sizeA directly.
     *
     * @param tGT   gene tree
     * @param tC    exemplar tree for cluster A
     * @param loC, hiC  range of A in tC's postorder array (the un-complemented range)
     * @param cComp true if A is a complement cluster (A = S\[loC,hiC))
     */
    public static int intersectWithFullTree(Tree tGT, Tree tC, int loC, int hiC, boolean cComp) {
        int L_GT = tGT.leafCount;
        int core = coreIntersect(tGT, 0, L_GT, tC, loC, hiC);
        return cComp ? (L_GT - core) : core;
    }

    // -------------------------------------------------------------------------
    // Multi-range cluster variants.
    //
    // A multi-range cluster's positive part is a union of PAIRWISE-DISJOINT
    // ranges {[los[j],his[j])} in its exemplar tree tC.  Because the ranges are
    // disjoint, |M ∩ (⋃ ranges)| = Σ_j |M ∩ [los[j],his[j])| — so each core count
    // is a sum of per-range coreIntersect()s (each still walks the smaller side).
    // The complement subtract trick (cComp ? size−core : core) is unchanged:
    // |comp ∩ M| = |M| − |(⋃ ranges) ∩ M|, valid since M ⊆ Lg ⊆ S.
    //
    // These mirror the single-range methods exactly when los.length == 1, so a
    // caller may always dispatch on Cluster.isMultiRange() with identical results.
    // -------------------------------------------------------------------------

    /** Σ_j |[loGT,hiGT) ∩ [los[j],his[j])| — disjoint ranges, each smaller-side walk. */
    public static int coreIntersectMulti(Tree tGT, int loGT, int hiGT,
                                          Tree tC, int[] los, int[] his) {
        int total = 0;
        for (int j = 0; j < los.length; j++) {
            total += coreIntersect(tGT, loGT, hiGT, tC, los[j], his[j]);
        }
        return total;
    }

    /** Multi-range analogue of {@link #intersect}. */
    public static int intersectMulti(Tree tGT, int loGT, int hiGT,
                                     Tree tC, int[] los, int[] his, boolean cComp,
                                     int sizeGTRange) {
        int core = coreIntersectMulti(tGT, loGT, hiGT, tC, los, his);
        return cComp ? (sizeGTRange - core) : core;
    }

    /** Multi-range analogue of {@link #intersectWithFullTree} (row sum vs full gene tree). */
    public static int intersectWithFullTreeMulti(Tree tGT, Tree tC,
                                                 int[] los, int[] his, boolean cComp) {
        int L_GT = tGT.leafCount;
        int core = coreIntersectMulti(tGT, 0, L_GT, tC, los, his);
        return cComp ? (L_GT - core) : core;
    }
}
