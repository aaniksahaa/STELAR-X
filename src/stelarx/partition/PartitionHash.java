package stelarx.partition;

import stelarx.cluster.ClusterHash;

/**
 * Order-invariant hash key for a gene-tree d-partition.
 *
 * BINARY (d=3): a rooted child partition is identified by the unordered pair
 * {M1,M2}. Taxa outside the node are deliberately absent from the key: the
 * rooted-triplet objective depends only on the actual children.
 * This matters for incomplete trees, where the same child partition can occur under
 * different ambient leaf sets.
 *
 * POLYTOMOUS (d≥4): a rooted node is identified by the unordered multiset of
 * its child hashes. No complement is accepted or stored.
 *
 * The two representations never collide: binary partitions (d=3) and polytomous ones
 * (d≥4) come from disjoint parser paths and {@code equals} short-circuits on {@code d}.
 */
public final class PartitionHash {

    private final int cachedHashCode;
    private final int d;       // number of parts (3 for binary, k+1 for polytomous)

    /** Flattened canonical fingerprints for both binary and general partitions. */
    private final long[] data;

    public PartitionHash(ClusterHash a, ClusterHash b) {
        // Decide ordering of the two explicit parts (a,b are interchangeable)
        boolean aFirst = compare(a, b) <= 0;
        ClusterHash first  = aFirst ? a : b;
        ClusterHash second = aFirst ? b : a;

        int m = a.sums.length;
        this.d = 3;
        data = new long[4 * m + 2];
        for (int s = 0; s < m; s++) {
            data[s]         = first.sums[s];
            data[s + m]     = first.xors[s];
            data[s + 2 * m] = second.sums[s];
            data[s + 3 * m] = second.xors[s];
        }
        data[4 * m] = first.size;
        data[4 * m + 1] = second.size;

        int h = 1;
        for (long v : data) h = 31 * h + Long.hashCode(v);
        this.cachedHashCode = h;
    }

    /**
     * General rooted-polytomy key, order-invariant over all supplied children.
     */
    public PartitionHash(ClusterHash[] children) {
        int childCount = children.length;
        int m = children[0].sums.length;
        long[][] fps = new long[childCount][2 * m + 1];
        for (int i = 0; i < childCount; i++) {
            for (int s = 0; s < m; s++) {
                fps[i][s]     = children[i].sums[s];
                fps[i][s + m] = children[i].xors[s];
            }
            fps[i][2 * m] = children[i].size;
        }
        java.util.Arrays.sort(fps, (x, y) -> {
            for (int s = 0; s < x.length; s++) {
                int c = Long.compareUnsigned(x[s], y[s]);
                if (c != 0) return c;
            }
            return 0;
        });
        long[] flat = new long[childCount * (2 * m + 1)];
        int p = 0;
        for (int i = 0; i < childCount; i++)
            for (long v : fps[i]) flat[p++] = v;

        this.d = childCount + 1;
        this.data = flat;
        int h = 1;
        for (long v : flat) h = 31 * h + Long.hashCode(v);
        this.cachedHashCode = h;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PartitionHash p)) return false;
        if (d != p.d) return false;
        return java.util.Arrays.equals(data, p.data);
    }

    @Override
    public int hashCode() { return cachedHashCode; }

    /** Lexicographic comparison of two ClusterHash objects (sums first, then xors). */
    private static int compare(ClusterHash a, ClusterHash b) {
        int m = a.sums.length;
        for (int s = 0; s < m; s++) {
            int c = Long.compareUnsigned(a.sums[s], b.sums[s]);
            if (c != 0) return c;
        }
        for (int s = 0; s < m; s++) {
            int c = Long.compareUnsigned(a.xors[s], b.xors[s]);
            if (c != 0) return c;
        }
        return Integer.compare(a.size, b.size);
    }
}
