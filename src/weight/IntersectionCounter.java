package weight;

import tree.Tree;

/**
 * O(min(|P|, |Q|)) intersection counts between postorder ranges.
 *
 * For complete gene trees (Lg = L), complement intersections reduce to:
 *   |comp(P) n Q| = |Q| - |P n Q|
 *   |P n comp(Q)| = |P| - |P n Q|
 *   |comp(P) n comp(Q)| = n - |P| - |Q| + |P n Q|
 *
 * So we only ever need the one "core" non-complement count.
 */
public final class IntersectionCounter {

    private IntersectionCounter() {}

    /**
     * |range_in_treeA n range_in_treeB| -- both ranges non-complement.
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
     * |M_range n cluster| where M_range is non-complement (a gene-tree subtree range)
     * and cluster may be complement (a species-tree candidate part).
     *
     * For complete trees: |comp(cluster_range) n M_range| = |M_range| - |cluster_range n M_range|
     *
     * @param tGT        gene-tree tree
     * @param loGT, hiGT range in gene-tree postorder array (non-complement)
     * @param tC         exemplar tree of the candidate cluster
     * @param loC, hiC   range in candidate cluster's exemplar tree
     * @param cComp      whether the candidate cluster is complement w.r.t. its tree
     * @param sizeGTRange = hiGT - loGT (passed in to avoid recomputation)
     */
    public static int intersect(Tree tGT, int loGT, int hiGT,
                                 Tree tC, int loC, int hiC, boolean cComp,
                                 int sizeGTRange) {
        int core = coreIntersect(tGT, loGT, hiGT, tC, loC, hiC);
        return cComp ? (sizeGTRange - core) : core;
    }

    /**
     * |A n Lg_GT| -- how many of the gene-tree's leaves fall inside cluster A.
     *
     * For a super-complement cluster A = S\sub(u):
     *   |A n Lg_GT| = L_GT - |sub(u) n Lg_GT|
     * For a plain cluster A = sub(u):
     *   |A n Lg_GT| = |sub(u) n Lg_GT|
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
}
