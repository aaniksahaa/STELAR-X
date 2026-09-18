import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SimilarityMatrixBuilder;
import stelarx.completion.SortedRowsBuilder;
import stelarx.completion.TreeCompleter;
import stelarx.completion.UPGMAClusterer;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.List;

/** Forces the large-N representation on small inputs and checks exact parity. */
public final class PackedSimilarityParityTest {
    private static void equalBits(double expected, double actual, String label) {
        if (Double.doubleToRawLongBits(expected) != Double.doubleToRawLongBits(actual)) {
            throw new AssertionError(label + ": expected=" + expected + " actual=" + actual);
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass one incomplete tree file");
        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        try {
            run(args[0]);
        } finally {
            Threading.shutdown();
        }
    }

    private static void run(String path) throws Exception {
        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(path, registry);
        int n = registry.size();

        System.clearProperty("stelarx.similarity.forcePacked");
        SimilarityMatrix dense = SimilarityMatrixBuilder.buildCPU(trees, n);
        System.setProperty("stelarx.similarity.forcePacked", "true");
        SimilarityMatrix packed;
        try {
            packed = SimilarityMatrixBuilder.buildCPU(trees, n);
        } finally {
            System.clearProperty("stelarx.similarity.forcePacked");
        }
        if (!packed.isPacked() || dense.isPacked()) throw new AssertionError("storage dispatch failed");

        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                equalBits(dense.getSim(a, b), packed.getSim(a, b), "similarity " + a + "," + b);
                equalBits(dense.getDist(a, b), packed.getDist(a, b), "distance " + a + "," + b);
            }
        }

        int[] denseRows = SortedRowsBuilder.buildCPU(dense.dist, n);
        int[][] packedRows = SortedRowsBuilder.buildPackedCPU(packed);
        for (int a = 0; a < n; a++) {
            for (int rank = 0; rank < n; rank++) {
                if (denseRows[a * n + rank] != packedRows[a][rank]) {
                    throw new AssertionError("sorted-row mismatch at row=" + a + " rank=" + rank);
                }
            }
        }

        Tree denseGuide = UPGMAClusterer.build(dense, trees.size());
        Tree packedGuide = UPGMAClusterer.build(packed, trees.size());
        String denseGuideNewick = denseGuide.toNewick(registry);
        String packedGuideNewick = packedGuide.toNewick(registry);
        if (!denseGuideNewick.equals(packedGuideNewick)) {
            throw new AssertionError("UPGMA mismatch:\n" + denseGuideNewick + "\n" + packedGuideNewick);
        }

        List<Tree> denseCompleted = TreeCompleter.completeAll(trees, dense, n);
        List<Tree> packedCompleted = TreeCompleter.completeAll(trees, packed, n);
        for (int i = 0; i < trees.size(); i++) {
            String expected = denseCompleted.get(i).toNewick(registry);
            String actual = packedCompleted.get(i).toNewick(registry);
            if (!expected.equals(actual)) {
                throw new AssertionError("completion mismatch for tree " + i
                    + ":\n" + expected + "\n" + actual);
            }
        }
        System.out.println("Packed similarity/storage/UPGMA/completion parity: PASS");
    }
}
