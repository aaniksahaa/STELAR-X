package cluster;

/**
 * Integer-tuple representation of a cluster (set of taxa).
 *
 * A cluster is always a range or complement-of-range in some gene tree's
 * postorder array, so 4 integers capture it exactly with O(1) memory.
 *
 *   complement=false:  taxa in postorderArray[treeIdx][left .. right)
 *   complement=true:   taxa in ALL taxa set NOT in [left..right)
 *                      i.e. super-complement: S \ [left..right)
 *
 * Size is pre-computed: either (right-left) or (totalTaxa - (right-left)).
 */
public final class Cluster {
    public final int treeIndex;
    public final int left;        // inclusive range start
    public final int right;       // exclusive range end
    public final boolean complement;
    public final int size;        // number of taxa in this cluster

    public Cluster(int treeIndex, int left, int right, boolean complement, int leafCount) {
        this.treeIndex  = treeIndex;
        this.left       = left;
        this.right      = right;
        this.complement = complement;
        int rangeSize   = right - left;
        this.size       = complement ? (leafCount - rangeSize) : rangeSize;
    }

    @Override
    public String toString() {
        return String.format("Cluster{t=%d, [%d,%d)%s, sz=%d}",
            treeIndex, left, right, complement ? "^c" : "", size);
    }
}
