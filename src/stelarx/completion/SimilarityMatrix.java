package stelarx.completion;

/**
 * Quartet-based taxon similarity matrix.
 *
 * For taxa a and b:
 *   M[a][b] = Σ_{t: a,b present} num_t(a,b)
 *             ──────────────────────────────
 *             Σ_{t: a,b present} den_t(a,b)
 *
 * where:
 *   num_t(a,b) = S[LCA_t(a,b)] − C2(subLC[c_a]) − C2(subLC[c_b])
 *   den_t(a,b) = C2(k_t − 2)
 *   k_t        = leaf count of tree t
 *
 * After normalize():
 *   sim[a][b]  ∈ [0,1], diagonal = 1, symmetric
 *   dist[a][b] = 1 − sim[a][b]   (used by TreeCompleter)
 *
 * Pairs that never co-occur in any tree have sim = 0, dist = 1.
 */
public class SimilarityMatrix {
    private static final int MAX_JAVA_ARRAY_LENGTH = Integer.MAX_VALUE - 8;

    public final int n;

    /** True when the exact symmetric matrix is held as segmented upper triangles. */
    private final boolean packed;

    /** Accumulated numerator: Σ num_t(a,b). Flat n×n double. */
    final double[] numSum;

    /** Accumulated denominator: Σ den_t(a,b). Flat n×n double. */
    final double[] denSum;

    /** Finalized similarity: numSum / denSum. Populated by normalize(). */
    public final double[] sim;

    /**
     * Finalized distance: 1 − sim.
     * Passed directly to TreeCompleter.completeAll() as the dist array.
     */
    public final double[] dist;

    private SegmentedDoubleArray packedNum;
    private SegmentedDoubleArray packedDen;
    /** True only when similarity preprocessing crossed the wide Java-array boundary. */
    private boolean streamedHostBatches;

    public SimilarityMatrix(int n) {
        long cellsLong = (long)n * n;
        this.n = n;
        this.packed = Boolean.getBoolean("stelarx.similarity.forcePacked")
            || requiresPacked(n);
        if (packed) {
            long triangleCells = triangleCellCount(n);
            this.numSum = null;
            this.denSum = null;
            this.sim = null;
            this.dist = null;
            this.packedNum = new SegmentedDoubleArray(triangleCells);
            this.packedDen = new SegmentedDoubleArray(triangleCells);
        } else {
            int cells = (int)cellsLong;
            this.numSum = new double[cells];
            this.denSum = new double[cells];
            this.sim = new double[cells];
            this.dist = new double[cells];
        }
    }

    /**
     * Finalize: compute sim[a][b] = numSum / denSum, then dist = 1 − sim.
     * Pairs that never co-occur get sim = 0, dist = 1.
     * Diagonal is set to sim = 1, dist = 0.
     */
    public void normalize() {
        if (packed) {
            double[][] nums = packedNum.segments();
            double[][] dens = packedDen.segments();
            for (int s = 0; s < nums.length; s++) {
                double[] num = nums[s];
                double[] den = dens[s];
                for (int i = 0; i < num.length; i++) {
                    num[i] = den[i] > 0.0 ? num[i] / den[i] : 0.0;
                }
            }
            for (int a = 0; a < n; a++) packedNum.set(index(a, a), 1.0);
            // Make the (potentially multi-GiB) denominator immediately collectible.
            packedDen = null;
            return;
        }
        for (int i = 0; i < sim.length; i++) {
            sim [i] = (denSum[i] > 0.0) ? numSum[i] / denSum[i] : 0.0;
            dist[i] = 1.0 - sim[i];
        }
        // Diagonal: self-similarity = 1, self-distance = 0
        for (int a = 0; a < n; a++) {
            sim [a * n + a] = 1.0;
            dist[a * n + a] = 0.0;
        }
    }

    /** Returns M[a][b]. Call after normalize(). */
    public double getSim(int a, int b)  {
        return packed ? packedNum.get(index(a, b)) : sim[a * n + b];
    }

    /** Returns (1 − M[a][b]). Call after normalize(). */
    public double getDist(int a, int b) {
        return packed ? 1.0 - packedNum.get(index(a, b)) : dist[a * n + b];
    }

    public boolean isPacked() { return packed; }

    void markStreamedHostBatches() { streamedHostBatches = true; }
    public boolean usedStreamedHostBatches() { return streamedHostBatches; }

    static long triangleCellCount(int n) {
        return (long)n * (n + 1L) / 2L;
    }

    static boolean requiresPacked(int n) {
        return (long)n * n > MAX_JAVA_ARRAY_LENGTH;
    }

    static long packedIndex(int n, int a, int b) {
        if (a > b) { int t = a; a = b; b = t; }
        return (long)a * n - (long)a * (a + 1L) / 2L + b;
    }

    long index(int a, int b) {
        return packedIndex(n, a, b);
    }

    void addPackedNumerator(int a, int b, double value) {
        packedNum.add(index(a, b), value);
    }

    void addPackedDenominator(int a, int b, double value) {
        packedDen.add(index(a, b), value);
    }

    double[][] packedNumeratorSegments() { return packedNum.segments(); }
    double[][] packedDenominatorSegments() { return packedDen.segments(); }
    int packedSegmentShift() { return SegmentedDoubleArray.SEGMENT_SHIFT; }

    SegmentedDoubleArray copyPackedSimilarity() {
        if (!packed || packedDen != null) {
            throw new IllegalStateException("packed similarity is not normalized");
        }
        SegmentedDoubleArray copy = new SegmentedDoubleArray(packedNum.length());
        double[][] src = packedNum.segments();
        double[][] dst = copy.segments();
        for (int s = 0; s < src.length; s++) {
            System.arraycopy(src[s], 0, dst[s], 0, src[s].length);
        }
        return copy;
    }
}
