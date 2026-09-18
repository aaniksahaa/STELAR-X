package stelarx.completion;

import stelarx.Config;

/** Arithmetic-only checks at and beyond the Java dense-array boundary. */
public final class PackedMatrixBoundaryTest {
    public static void main(String[] args) {
        if (SimilarityMatrix.requiresPacked(46_340)) {
            throw new AssertionError("46,340 should retain the established dense path");
        }
        if (!SimilarityMatrix.requiresPacked(46_341)
                || !SimilarityMatrix.requiresPacked(50_000)) {
            throw new AssertionError("large-N packed dispatch boundary is wrong");
        }

        int n = 50_000;
        long cells = SimilarityMatrix.triangleCellCount(n);
        if (cells != 1_250_025_000L) throw new AssertionError("triangle size=" + cells);
        if (SimilarityMatrix.packedIndex(n, 0, 0) != 0L
                || SimilarityMatrix.packedIndex(n, 0, n - 1) != n - 1L
                || SimilarityMatrix.packedIndex(n, n - 1, n - 1) != cells - 1L) {
            throw new AssertionError("packed boundary indices are wrong");
        }
        int[][] pairs = {{1, 49_999}, {23_456, 45_678}, {49_998, 49_999}};
        for (int[] p : pairs) {
            long ab = SimilarityMatrix.packedIndex(n, p[0], p[1]);
            long ba = SimilarityMatrix.packedIndex(n, p[1], p[0]);
            if (ab != ba || ab < 0 || ab >= cells) {
                throw new AssertionError("invalid symmetric index for " + p[0] + "," + p[1]);
            }
        }

        SegmentedDoubleArray segmented = new SegmentedDoubleArray(40, 4);
        if (segmented.segments().length != 3
                || segmented.segments()[0].length != 16
                || segmented.segments()[1].length != 16
                || segmented.segments()[2].length != 8) {
            throw new AssertionError("test segment layout is wrong");
        }
        long[] positions = {0, 15, 16, 31, 32, 39};
        for (long p : positions) segmented.set(p, p + 0.25);
        segmented.add(16, 2.0);
        for (long p : positions) {
            double expected = p + 0.25 + (p == 16 ? 2.0 : 0.0);
            if (segmented.get(p) != expected) {
                throw new AssertionError("segment-boundary access failed at " + p);
            }
        }

        System.setProperty("stelarx.similarity.forcePacked", "true");
        SimilarityMatrix tinyPacked;
        try {
            tinyPacked = new SimilarityMatrix(3);
        } finally {
            System.clearProperty("stelarx.similarity.forcePacked");
        }
        Config cfg = Config.getInstance();
        if (SimilarityMatrixBuilder.effectiveTreeCapMiB(tinyPacked, cfg) != 8192) {
            throw new AssertionError("large-N automatic GPU batching ceiling was not raised");
        }
        cfg.setGpuSimilarityVramCapMiB(256);
        if (SimilarityMatrixBuilder.effectiveTreeCapMiB(tinyPacked, cfg) != 256) {
            throw new AssertionError("explicit GPU batching ceiling was not respected");
        }

        // For the reported 1,000-taxon trees E=2,998. The former 75k case must
        // retain its established one-shot wide path, while 100k must stream
        // because its wide micro-RMQ array would exceed Java's element limit.
        long arrayLimit = Integer.MAX_VALUE - 8L;
        if (!SimilarityMatrixBuilder.wideLayoutFits(75_000, 1_000, 2_998, arrayLimit)) {
            throw new AssertionError("75k wide layout was unnecessarily moved off one-shot mode");
        }
        if (SimilarityMatrixBuilder.wideLayoutFits(100_000, 1_000, 2_998, arrayLimit)) {
            throw new AssertionError("100k overflowing wide layout was accepted as one-shot");
        }

        int streamed = SimilarityMatrixBuilder.wideBatchTreeCount(
            100_000, 1_000, 2_998, arrayLimit, 1L << 30);
        if (streamed != 7_718) {
            throw new AssertionError("100k wide batch size changed: " + streamed);
        }
        long microCells = (long)streamed * 9 * 2_998;
        long flatBytes = (long)streamed
            * (36L * 2_998 + 9L * 2_998 + 4L * 4 * 12 + 4L * 1_000 + 8L);
        if (microCells > arrayLimit || flatBytes > (1L << 30)) {
            throw new AssertionError("streamed wide batch exceeds a safety bound");
        }
        long nextFlatBytes = (long)(streamed + 1)
            * (36L * 2_998 + 9L * 2_998 + 4L * 4 * 12 + 4L * 1_000 + 8L);
        if (nextFlatBytes <= (1L << 30)) {
            throw new AssertionError("streamed wide batch is smaller than necessary");
        }

        // Small fitting inputs retain one batch, and artificial tiny limits
        // exercise the exact per-array planner boundary without large arrays.
        if (SimilarityMatrixBuilder.wideBatchTreeCount(
                7, 10, 16, 2_000, 1_000_000) != 7) {
            throw new AssertionError("fitting wide input was unnecessarily split");
        }
        if (SimilarityMatrixBuilder.wideBatchTreeCount(
                100, 10, 16, 1_000, 1_000_000) != 6) {
            throw new AssertionError("wide host-byte boundary was not enforced exactly");
        }
        System.out.println("Packed matrix 46,340/46,341/50,000 boundaries: PASS");
    }
}
