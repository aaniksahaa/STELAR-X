package partition;

import cluster.ClusterHash;

/**
 * A gene-tree tripartition for STELAR-X triplet scoring.
 *
 * For a non-root internal node u of a rooted gene tree g:
 *   part1 = left-subtree  of u  (range [leftStart, leftEnd))
 *   part2 = right-subtree of u  (range [rightStart, rightEnd))
 *   part3 = Lg \ (part1 u part2)  -- "everything else in this tree"
 *
 * Sizes of all three parts are stored for O(1) triplet score computation.
 */
public final class RangePartition {

    /** Hash of part1 (left subtree). */
    public final ClusterHash hash1;

    /** Hash of part2 (right subtree). */
    public final ClusterHash hash2;

    /** Hash of part3 = Lg \ (part1 u part2). */
    public final ClusterHash hash3;

    /** Sizes of the three parts. */
    public final int size1, size2, size3;

    /** Index of the gene tree this partition came from. */
    public final int treeIndex;

    /** Range of the left child in the exemplar tree's postorderArray. */
    public final int leftStart, leftEnd;

    /** Range of the right child in the exemplar tree's postorderArray. */
    public final int rightStart, rightEnd;

    public RangePartition(ClusterHash h1, ClusterHash h2, ClusterHash h3,
                          int sz1, int sz2, int sz3,
                          int treeIndex,
                          int leftStart, int leftEnd,
                          int rightStart, int rightEnd) {
        this.hash1 = h1; this.hash2 = h2; this.hash3 = h3;
        this.size1 = sz1; this.size2 = sz2; this.size3 = sz3;
        this.treeIndex = treeIndex;
        this.leftStart = leftStart; this.leftEnd = leftEnd;
        this.rightStart = rightStart; this.rightEnd = rightEnd;
    }

    @Override
    public String toString() {
        return String.format("RangePartition{sz=%d|%d|%d, t=%d, L=[%d,%d) R=[%d,%d)}",
            size1, size2, size3, treeIndex,
            leftStart, leftEnd, rightStart, rightEnd);
    }
}
