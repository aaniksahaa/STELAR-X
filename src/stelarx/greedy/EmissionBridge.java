package stelarx.greedy;

import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.tree.Tree;

import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.List;

/**
 * Bridges polytomy-resolution emissions into the DP cluster set X.
 *
 * Each {@link EmittedBipartition} carries a {@code ClusterHash} signature (computed
 * via the consensus tree's prefix scans with the SAME {@code TaxonHasher} as the
 * gene-tree clusters, so it is directly comparable) and a {@link MultiRange}
 * descriptor of its smaller side in some consensus snapshot's {@code aCons}.
 *
 * Two tiers:
 *   Tier 1 — the set is already in X (some gene tree realizes it contiguously):
 *            the signature is present, so nothing to synthesize — it is already a
 *            scorable single-range cluster.
 *   Tier 2 — the set is contiguous in no tree: synthesize a multi-range
 *            {@link Cluster} anchored on the consensus snapshot tree (a COMPLETE
 *            exemplar over all n taxa) and insert it into X.
 *
 * Consensus snapshot trees that contribute a synthesized exemplar are wrapped as
 * lightweight {@link Tree}s (postorderArray = {@code aCons}, inverse positionMap,
 * {@code root == null} since they are membership-only) and APPENDED to a
 * weight-calc exemplar list — NOT to the gene-tree list the DP mines for local
 * transitions. The synthesized cluster's {@code treeIndex} points into that
 * appended region. Mode 2 (cross-tree) then discovers the transitions that make
 * the cluster usable in the DP (it is hash-based, so it incorporates the new
 * clusters automatically).
 */
public final class EmissionBridge {

    private EmissionBridge() {}

    /**
     * Insert all emissions into {@code clusterTable}, appending any needed
     * consensus exemplar trees to {@code clusterTreesMutable}.
     *
     * @param clusterTreesMutable  the WEIGHT-calc exemplar list (a copy of the
     *                             gene/completed-tree list); this method appends
     *                             consensus exemplar trees to it. The DP's own
     *                             gene-tree list must NOT be this list.
     * @return {@code int[]{tier1, tier2}} — counts of already-present vs synthesized.
     */
    public static int[] bridge(EmissionBuffer buffer, ClusterTable clusterTable,
                               List<Tree> clusterTreesMutable, int numTaxa) {
        return bridge(buffer, clusterTable, clusterTreesMutable, numTaxa, false, -1);
    }

    /**
     * @param anchorFreeX  when true (anchor-free X), register each emission's
     *                     ANCHOR-FREE orientation — the side not containing the anchor.
     *                     The canonical side {@code E} is registered as-is when it is
     *                     anchor-free, else its complement {@code S\E} is synthesized
     *                     (multi-range, complement=true; the consensus exemplar spans all
     *                     n taxa, so the complement is walkable). Keeps X orientation-
     *                     complete under anchoring so the DP can still use the emission.
     * @param anchor       anchor taxon global id (valid iff anchorFreeX).
     */
    public static int[] bridge(EmissionBuffer buffer, ClusterTable clusterTable,
                               List<Tree> clusterTreesMutable, int numTaxa,
                               boolean anchorFreeX, int anchor) {
        int tier1 = 0, tier2 = 0;
        // One exemplar Tree per distinct consensus snapshot that needs synthesis.
        IdentityHashMap<ConsensusTree, Integer> exemplarIndex = new IdentityHashMap<>();
        // Anchor position per consensus tree (cached; consensus trees span all n taxa).
        IdentityHashMap<ConsensusTree, Integer> anchorPosCache =
            anchorFreeX ? new IdentityHashMap<>() : null;
        ClusterHash root = clusterTable.getAllTaxaHash();

        for (EmittedBipartition e : buffer.all()) {
            MultiRange mr = e.canonicalSide;
            ConsensusTree ct = mr.tree;

            // Decide which orientation to register.  Default: the canonical side E.
            ClusterHash regHash   = e.signature;
            boolean     regCompl  = false;
            int         regSize   = e.size;
            if (anchorFreeX) {
                int aPos = anchorPosCache.computeIfAbsent(ct, k -> anchorPosIn(k, anchor));
                if (inRanges(aPos, mr.los, mr.his)) {   // E contains the anchor → use S\E
                    regHash  = ClusterHash.residual(root, e.signature);
                    regCompl = true;
                    regSize  = numTaxa - e.size;
                }
            }

            // Tier 1: the orientation we want is already a scorable cluster.
            if (clusterTable.contains(regHash)) { tier1++; continue; }

            // Tier 2: synthesize a consensus-anchored multi-range exemplar.
            Integer ti = exemplarIndex.get(ct);
            if (ti == null) {
                ti = clusterTreesMutable.size();
                clusterTreesMutable.add(buildExemplar(ct, ti, numTaxa));
                exemplarIndex.put(ct, ti);
            }
            Cluster c = new Cluster(ti, mr.los, mr.his, regCompl, regSize);
            if (clusterTable.addCluster(regHash, c)) tier2++;
        }
        return new int[]{ tier1, tier2 };
    }

    /** Postorder index of the anchor taxon in a consensus tree (spans all n taxa). */
    private static int anchorPosIn(ConsensusTree ct, int anchor) {
        int[] aCons = ct.aCons();
        for (int p = 0; p < aCons.length; p++) if (aCons[p] == anchor) return p;
        return -1;
    }

    /** True iff pos falls in any half-open range [los[j],his[j]). */
    private static boolean inRanges(int pos, int[] los, int[] his) {
        if (pos < 0) return false;
        for (int j = 0; j < los.length; j++) if (pos >= los[j] && pos < his[j]) return true;
        return false;
    }

    /**
     * Wrap a consensus snapshot as a membership-only {@link Tree}: postorderArray
     * = {@code aCons}, inverse positionMap, {@code root == null}. The consensus
     * tree spans all n taxa, so {@code isComplete} is true and complement-free
     * intersection works directly.
     */
    private static Tree buildExemplar(ConsensusTree ct, int treeIndex, int numTaxa) {
        int[] aCons = ct.aCons();
        int[] posMap = new int[numTaxa];
        Arrays.fill(posMap, -1);
        for (int p = 0; p < aCons.length; p++) posMap[aCons[p]] = p;
        return new Tree(treeIndex, /*root=*/null, aCons, posMap, aCons.length, numTaxa);
    }
}
