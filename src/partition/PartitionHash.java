package partition;

import cluster.ClusterHash;

/**
 * Order-invariant hash key for a gene-tree tripartition.
 *
 * A tripartition (M1|M2|M3) is identified by the unordered pair {M1,M2}
 * plus M3. The combined hash is order-invariant over (h1,h2) and then includes h3.
 */
public final class PartitionHash {

    private final int cachedHashCode;

    private final long[] lo;   // the lexicographically smaller of {h1,h2}
    private final long[] hi;   // the larger
    private final long[] m3;   // h3 (M3 = Lg \ M1 \ M2)

    public PartitionHash(ClusterHash a, ClusterHash b, ClusterHash c) {
        boolean aFirst = compare(a, b) <= 0;
        ClusterHash first  = aFirst ? a : b;
        ClusterHash second = aFirst ? b : a;

        int m = a.sums.length;
        lo = new long[2 * m];
        hi = new long[2 * m];
        m3 = new long[2 * m];
        for (int s = 0; s < m; s++) {
            lo[s]     = first.sums[s];
            lo[s + m] = first.xors[s];
            hi[s]     = second.sums[s];
            hi[s + m] = second.xors[s];
            m3[s]     = c.sums[s];
            m3[s + m] = c.xors[s];
        }

        int h = 1;
        for (long v : lo) h = 31 * h + Long.hashCode(v);
        for (long v : hi) h = 31 * h + Long.hashCode(v);
        for (long v : m3) h = 31 * h + Long.hashCode(v);
        this.cachedHashCode = h;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PartitionHash)) return false;
        PartitionHash p = (PartitionHash) o;
        if (lo.length != p.lo.length) return false;
        for (int i = 0; i < lo.length; i++) {
            if (lo[i] != p.lo[i] || hi[i] != p.hi[i] || m3[i] != p.m3[i]) return false;
        }
        return true;
    }

    @Override
    public int hashCode() { return cachedHashCode; }

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
        return 0;
    }
}
