package stelarx.greedy;

/**
 * One side of a polytomy-derived bipartition, represented as a union of
 * pairwise-disjoint half-open ranges {@code [los[j], his[j])} into the parent
 * {@link ConsensusTree}'s {@code aCons} array.
 *
 * The disjointness precondition matters: ϕ1 (additive) double-counts overlaps
 * and ϕ2 (XOR) silently cancels duplicated elements (design §7.3).  All
 * emissions in the polytomy-resolution phase combine consensus-tree CHILD
 * ranges of a single polytomy node, which are disjoint by construction —
 * but new callers MUST preserve that invariant.
 */
public final class MultiRange {
    public final ConsensusTree tree;
    public final int[] los;
    public final int[] his;

    public MultiRange(ConsensusTree tree, int[] los, int[] his) {
        if (los.length != his.length)
            throw new IllegalArgumentException("MultiRange los/his length mismatch");
        this.tree = tree;
        this.los  = los;
        this.his  = his;
    }

    /** Total number of taxa across all ranges. */
    public int size() {
        int sz = 0;
        for (int j = 0; j < los.length; j++) sz += his[j] - los[j];
        return sz;
    }

    public int numRanges() { return los.length; }
}
