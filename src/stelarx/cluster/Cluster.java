package stelarx.cluster;

/**
 * Integer-tuple representation of a cluster (set of taxa).
 *
 * <h3>Single-range (the common, fast case)</h3>
 * A cluster is a range or complement-of-range in some tree's postorder array,
 * so 4 integers capture it exactly with O(1) memory.  {@code los/his == null}.
 *
 *   complement=false:  taxa in postorderArray[treeIdx][left .. right)
 *   complement=true:   taxa in that tree NOT in [left..right)
 *                      i.e. postorderArray[treeIdx][0..left) ∪ [right..leafCount)
 *
 * Size is pre-computed: either (right-left) or (leafCount - (right-left)).
 *
 * <h3>Multi-range (the fallback for non-contiguous sets)</h3>
 * Some sets — e.g. polytomy-resolution emissions (consensus-design §7.3) — are
 * contiguous in no tree.  They are stored as a union of pairwise-DISJOINT ranges
 * {@code [los[j], his[j])} in an exemplar tree's postorder ({@code los != null}).
 * The exemplar must be a COMPLETE tree (all taxa positioned), e.g. a consensus
 * tree, so the set and its complement are both expressible.
 * {@code complement=true} then means
 * S \ (⋃ ranges), handled by the same subtract-from-|M| identity as single-range.
 * {@code left/right} hold the bounding span [los[0], his[last]) for debug only;
 * never used for membership when {@code los != null}.
 */
public final class Cluster {
    public final int treeIndex;
    public final int left;        // inclusive range start (single-range) / span lo (multi)
    public final int right;       // exclusive range end   (single-range) / span hi (multi)
    public final boolean complement;
    public final int size;        // number of taxa in this cluster

    /** Multi-range descriptor: disjoint ranges [los[j],his[j]); null for single-range. */
    public final int[] los;
    public final int[] his;

    /** Single-range cluster (los/his == null) — the existing fast path, unchanged. */
    public Cluster(int treeIndex, int left, int right, boolean complement, int treeLeafCount) {
        this.treeIndex  = treeIndex;
        this.left       = left;
        this.right      = right;
        this.complement = complement;
        int rangeSize   = right - left;
        this.size       = complement ? (treeLeafCount - rangeSize) : rangeSize;
        this.los = null;
        this.his = null;
    }

    /**
     * Multi-range cluster: union of disjoint ranges in {@code treeIndex}'s postorder.
     * @param los/his  disjoint half-open ranges (callers should pre-merge/sort)
     * @param complement  if true, the set is S \ (⋃ ranges)
     * @param size  precomputed |cluster| (caller knows it; e.g. emission size)
     */
    public Cluster(int treeIndex, int[] los, int[] his, boolean complement, int size) {
        if (los == null || his == null || los.length != his.length || los.length == 0)
            throw new IllegalArgumentException("multi-range Cluster needs equal non-empty los/his");
        this.treeIndex  = treeIndex;
        this.los        = los;
        this.his        = his;
        this.left       = los[0];
        this.right      = his[his.length - 1];
        this.complement = complement;
        this.size       = size;
    }

    /** True iff this cluster is stored as a union of ≥1 explicit ranges (los/his). */
    public boolean isMultiRange() { return los != null; }

    @Override
    public String toString() {
        if (los == null)
            return String.format("Cluster{t=%d, [%d,%d)%s, sz=%d}",
                treeIndex, left, right, complement ? "^c" : "", size);
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("Cluster{t=%d, multi[", treeIndex));
        for (int j = 0; j < los.length; j++) {
            if (j > 0) sb.append('∪');
            sb.append(String.format("[%d,%d)", los[j], his[j]));
        }
        sb.append(String.format("]%s, sz=%d}", complement ? "^c" : "", size));
        return sb.toString();
    }
}
