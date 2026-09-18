package stelarx.partition;

/**
 * A compact rooted gene-tree child partition.
 *
 * For an internal node with {@code k} actual children, only those child sets are
 * represented. Taxa elsewhere in the gene tree do not participate in the rooted
 * triplet objective and are deliberately absent from this data structure.
 *
 * {@code d = k + 1} is retained as a compact discriminator used throughout the
 * existing backends: {@code d == 3} means binary and {@code d > 3} means a native
 * polytomy. It must not be interpreted as saying that an outside set is stored.
 */
public final class Partition {

    /** Legacy discriminator: actual child count plus one. */
    public final int d;

    /** Actual child sizes and ranges; non-null only for polytomies, length d-1. */
    public final int[] sizes;
    public final int[] partStarts;
    public final int[] partEnds;

    /** Index of the gene tree this partition came from (of the first exemplar). */
    public final int treeIndex;

    // Binary (d=3) compact fields.
    public final int size1, size2;
    public final int leftStart, leftEnd, rightStart, rightEnd;

    /** Binary constructor. */
    public Partition(int sz1, int sz2, int treeIndex,
                     int leftStart, int leftEnd, int rightStart, int rightEnd) {
        this.d = 3;
        this.sizes = null;
        this.partStarts = null;
        this.partEnds = null;
        this.treeIndex = treeIndex;
        this.size1 = sz1;
        this.size2 = sz2;
        this.leftStart = leftStart;
        this.leftEnd = leftEnd;
        this.rightStart = rightStart;
        this.rightEnd = rightEnd;
    }

    /** Polytomy constructor; every array contains exactly the actual children. */
    public Partition(int[] sizes, int[] partStarts, int[] partEnds, int treeIndex) {
        if (sizes.length < 3 || partStarts.length != sizes.length || partEnds.length != sizes.length)
            throw new IllegalArgumentException("A polytomy requires matching arrays for at least 3 children");
        this.d = sizes.length + 1;
        this.sizes = sizes;
        this.partStarts = partStarts;
        this.partEnds = partEnds;
        this.treeIndex = treeIndex;
        this.size1 = sizes[0];
        this.size2 = sizes[1];
        this.leftStart = partStarts[0];
        this.leftEnd = partEnds[0];
        this.rightStart = partStarts[sizes.length - 1];
        this.rightEnd = partEnds[sizes.length - 1];
    }

    public boolean isPolytomous() { return d > 3; }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Partition{children=");
        if (d == 3) {
            sb.append(size1).append('|').append(size2);
        } else {
            for (int i = 0; i < sizes.length; i++) {
                if (i > 0) sb.append('|');
                sb.append(sizes[i]);
            }
        }
        return sb.append(", t=").append(treeIndex).append('}').toString();
    }
}
