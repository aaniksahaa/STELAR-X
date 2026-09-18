package stelarx.greedy;

import stelarx.cluster.ClusterHash;

/**
 * One bipartition emitted by polytomy resolution.
 *
 *   signature      — double-hash signature ({@code ClusterHash}) computed via
 *                    the consensus-tree prefix-scan; matches gene-tree
 *                    signatures for the same taxon set (design §7.2)
 *   canonicalSide  — the SMALLER side, as a multi-range descriptor pointing
 *                    into the source consensus tree's {@code aCons} array
 *   size           — |canonicalSide|
 *   source         — 'A' for Step A (UPGMA on group similarity), 'B' for Step B
 *                    (sampleAndResolve / induced-split sampling)
 *   thresholdIndex — which T[ti] snapshot the source polytomy lives in
 */
public final class EmittedBipartition {
    public final ClusterHash signature;
    public final MultiRange  canonicalSide;
    public final int         size;
    public final char        source;
    public final int         thresholdIndex;

    public EmittedBipartition(ClusterHash signature, MultiRange canonicalSide,
                              int size, char source, int thresholdIndex) {
        this.signature      = signature;
        this.canonicalSide  = canonicalSide;
        this.size           = size;
        this.source         = source;
        this.thresholdIndex = thresholdIndex;
    }
}
