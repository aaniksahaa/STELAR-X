package dp;

import cluster.ClusterHash;

/**
 * An order-invariant pair of ClusterHashes representing one candidate split.
 *
 * (a, b) and (b, a) are the same bipartition.
 * Canonical form: the lexicographically smaller hash is stored first.
 */
public final class BipartitionSplit {

    public final ClusterHash lo;   // canonical "smaller" half
    public final ClusterHash hi;   // canonical "larger"  half

    private final int cachedHashCode;

    public BipartitionSplit(ClusterHash a, ClusterHash b) {
        if (compare(a, b) <= 0) { this.lo = a; this.hi = b; }
        else                     { this.lo = b; this.hi = a; }
        this.cachedHashCode = 31 * lo.hashCode() + hi.hashCode();
    }

    // Lexicographic compare on (sums[], xors[]) arrays
    private static int compare(ClusterHash x, ClusterHash y) {
        for (int i = 0; i < x.sums.length; i++) {
            int c = Long.compareUnsigned(x.sums[i], y.sums[i]);
            if (c != 0) return c;
        }
        for (int i = 0; i < x.xors.length; i++) {
            int c = Long.compareUnsigned(x.xors[i], y.xors[i]);
            if (c != 0) return c;
        }
        return 0;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BipartitionSplit)) return false;
        BipartitionSplit s = (BipartitionSplit) o;
        return lo.equals(s.lo) && hi.equals(s.hi);
    }

    @Override
    public int hashCode() { return cachedHashCode; }

    @Override
    public String toString() {
        return "Split{" + lo + " | " + hi + "}";
    }
}
