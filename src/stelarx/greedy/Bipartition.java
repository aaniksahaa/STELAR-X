package stelarx.greedy;

import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;

/**
 * One unrooted bipartition selected from ClusterTable for greedy-consensus consumption.
 *
 * The bipartition has two sides; the {@code canonicalHash}/{@code canonicalExemplar}
 * fields point at the side chosen by {@link BipartitionCounter}'s canonicalization
 * rule.  The build phase enumerates this side's taxa to drive INSERT.
 *
 * {@code frequency} is the number of input gene trees in which this bipartition
 * appears (with either orientation).  Frequencies of the two sides in ClusterTable
 * are equal for completed-tree inputs, so we just use one.
 */
public final class Bipartition {
    public final ClusterHash canonicalHash;
    public final Cluster     canonicalExemplar;
    public final int         size;        // |canonical side|
    public final int         frequency;

    public Bipartition(ClusterHash canonicalHash, Cluster canonicalExemplar,
                       int size, int frequency) {
        this.canonicalHash     = canonicalHash;
        this.canonicalExemplar = canonicalExemplar;
        this.size              = size;
        this.frequency         = frequency;
    }
}
