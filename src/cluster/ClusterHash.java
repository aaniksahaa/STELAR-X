package cluster;

import java.util.Arrays;

/**
 * Immutable hash key for a cluster (set of taxa).
 *
 * Stores 2m values: m sum-hashes and m XOR-hashes.
 * Values are stored RAW (un-mixed) so that arithmetic operations are valid:
 *
 *   sum(A u B)  = sum(A)  + sum(B)   (mod 2^64)
 *   xor(A u B)  = xor(A)  ^ xor(B)   (for disjoint A, B)
 *   sum(A \ B)  = sum(A)  - sum(B)   (mod 2^64, valid when B subseteq A)
 *   xor(A \ B)  = xor(A)  ^ xor(B)   (valid when B subseteq A)
 *
 * Also stores the cluster size for O(1) access.
 */
public final class ClusterHash {

    public final long[] sums;   // sums[s]  = raw sum  of taxon hashes under seed s
    public final long[] xors;   // xors[s]  = raw XOR  of taxon hashes under seed s
    public final int size;      // number of taxa in this cluster
    private final int cachedHashCode;

    public ClusterHash(long[] rawSums, long[] rawXors, int size, int m) {
        this.size = size;
        this.sums = Arrays.copyOf(rawSums, m);
        this.xors = Arrays.copyOf(rawXors, m);
        int h = 1;
        for (long v : this.sums) h = 31 * h + Long.hashCode(v);
        for (long v : this.xors) h = 31 * h + Long.hashCode(v);
        this.cachedHashCode = h;
    }

    /**
     * Compute hash of set difference A \ B, given B subseteq A.
     * Result has sums[s] = a.sums[s] - b.sums[s] and xors[s] = a.xors[s] ^ b.xors[s].
     * Size = a.size - b.size.
     */
    public static ClusterHash residual(ClusterHash a, ClusterHash b) {
        int m = a.sums.length;
        long[] resSums = new long[m];
        long[] resXors = new long[m];
        for (int s = 0; s < m; s++) {
            resSums[s] = a.sums[s] - b.sums[s];
            resXors[s] = a.xors[s] ^ b.xors[s];
        }
        return new ClusterHash(resSums, resXors, a.size - b.size, m);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof ClusterHash)) return false;
        ClusterHash c = (ClusterHash) o;
        return size == c.size
            && Arrays.equals(sums, c.sums)
            && Arrays.equals(xors, c.xors);
    }

    @Override
    public int hashCode() { return cachedHashCode; }

    @Override
    public String toString() {
        return String.format("CH{size=%d, sum0=%016x, xor0=%016x}", size, sums[0], xors[0]);
    }
}
