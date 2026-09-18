package stelarx.completion;

/**
 * Average pairwise topological distance between taxa across gene trees.
 *
 * For taxa a and b, D[a][b] = average over all gene trees containing both a and b of:
 *   depth(a) + depth(b) - 2 * depth(LCA(a,b))
 *
 * This equals the number of edges on the path from a to b in that gene tree.
 *
 * Usage:
 *   1. Build via DistanceMatrixBuilder.buildCPU(trees, n).
 *   2. Query with get(a, b).
 *   3. Pairs that never co-occur in any tree return Double.MAX_VALUE.
 */
public class DistanceMatrix {
    public final int n;

    /** Accumulated sum of topological distances (before normalization). Flat n×n. */
    final double[] distSum;

    /** Number of gene trees where both taxa appear. Flat n×n. */
    final int[] cooccurrence;

    /** Finalized average distances. Populated by normalize(). Flat n×n. */
    public final double[] dist;

    public DistanceMatrix(int n) {
        this.n            = n;
        this.distSum      = new double[n * n];
        this.cooccurrence = new int[n * n];
        this.dist         = new double[n * n];
    }

    /**
     * Finalize: compute dist[a][b] = distSum[a][b] / cooccurrence[a][b].
     * Pairs that never co-occur receive Double.MAX_VALUE.
     * Self-distance is always 0.
     */
    public void normalize() {
        for (int i = 0; i < n * n; i++) {
            dist[i] = (cooccurrence[i] > 0)
                    ? distSum[i] / cooccurrence[i]
                    : Double.MAX_VALUE;
        }
        for (int a = 0; a < n; a++) dist[a * n + a] = 0.0;
    }

    /** Returns D[a][b]. Call after normalize(). */
    public double get(int a, int b) {
        return dist[a * n + b];
    }
}
