import stelarx.completion.EulerTourBuilder;
import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SimilarityMatrixBuilder;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;

import java.util.List;

/**
 * Exhaustively verifies the compact similarity-RMQ representation.
 *
 * For every power-of-two sparse interval and every arbitrary query interval,
 * the stored/queried position must equal a direct left-to-right scan's first
 * minimum-depth Euler position. Equality of the position is stronger than
 * equality of the minimum depth: it also proves that child-size/F payloads are
 * fetched from exactly the same left-biased Euler visit as before.
 */
public final class SimilarityArgminTest {
    private static boolean isStrictBinary(TreeNode node) {
        if (node.isLeaf()) return true;
        return !node.isPolytomous() && node.left != null && node.right != null
            && isStrictBinary(node.left) && isStrictBinary(node.right);
    }

    private static int directArgmin(short[] depths, int lo, int hiInclusive) {
        int best = lo;
        for (int i = lo + 1; i <= hiInclusive; i++) {
            if (depths[i] < depths[best]) best = i;
        }
        return best;
    }

    private static long verifyTree(Tree tree, int n) {
        EulerTourBuilder.FullTourData td = EulerTourBuilder.buildFull(tree, n);
        EulerTourBuilder.WideTourData wide = EulerTourBuilder.buildWide(tree, n);
        int len = td.tourLen;
        long checks = 0;

        if (wide.tourLen != len) {
            throw new AssertionError("wide tour length mismatch: compact=" + len
                + " wide=" + wide.tourLen);
        }
        for (int p = 0; p < len; p++) {
            if (wide.depths[p] != td.depths[p]
                    || wide.eulerF[p] != td.eulerF[p]
                    || wide.eulerLeftChildS[p] != td.eulerLeftChildS[p]
                    || wide.eulerLeftChildF[p] != td.eulerLeftChildF[p]
                    || wide.eulerRightChildS[p] != td.eulerRightChildS[p]
                    || wide.eulerRightChildF[p] != td.eulerRightChildF[p]) {
                throw new AssertionError("wide payload mismatch at Euler position " + p);
            }
        }
        for (int a = 0; a < n; a++) {
            if (wide.firstOcc[a] != td.firstOcc[a]) {
                throw new AssertionError("wide first-occurrence mismatch for taxon " + a);
            }
        }

        // Verify every materialized sparse-table cell independently.
        for (int lvl = 0; lvl < td.log; lvl++) {
            int width = 1 << lvl;
            int end = len - width + 1;
            for (int lo = 0; lo < end; lo++) {
                int expected = directArgmin(td.depths, lo, lo + width - 1);
                int actual = td.sparseArgmin[lvl][lo];
                if (actual != expected) {
                    throw new AssertionError("sparse cell mismatch: tree=" + tree.treeIndex
                        + " level=" + lvl + " lo=" + lo + " expected=" + expected
                        + " actual=" + actual);
                }
                checks++;
            }
        }

        // Verify the exact two-overlapping-block query used by the CUDA kernel.
        for (int lo = 0; lo < len; lo++) {
            for (int hi = lo; hi < len; hi++) {
                int width = hi - lo + 1;
                int lvl = 31 - Integer.numberOfLeadingZeros(width);
                int lo2 = hi - (1 << lvl) + 1;
                int posL = td.sparseArgmin[lvl][lo];
                int posR = td.sparseArgmin[lvl][lo2];
                int actual = (td.depths[posL] <= td.depths[posR]) ? posL : posR;
                int expected = directArgmin(td.depths, lo, hi);
                if (actual != expected) {
                    throw new AssertionError("query mismatch: tree=" + tree.treeIndex
                        + " interval=[" + lo + "," + hi + "] expected=" + expected
                        + " actual=" + actual);
                }
                int wideActual = EulerTourBuilder.queryWideArgmin(wide, lo, hi);
                if (wideActual != expected) {
                    throw new AssertionError("wide query mismatch: tree=" + tree.treeIndex
                        + " interval=[" + lo + "," + hi + "] expected=" + expected
                        + " actual=" + wideActual);
                }
                checks++;
            }
        }
        return checks;
    }

    /** Java emulation of the compact CUDA query; compare its final matrix to the CPU reference. */
    private static long verifyCompactSimilarity(List<Tree> trees, int n) {
        for (Tree tree : trees) if (!isStrictBinary(tree.root)) return 0;

        double[] num = new double[n * n];
        double[] den = new double[n * n];
        for (Tree tree : trees) {
            EulerTourBuilder.FullTourData td = EulerTourBuilder.buildFull(tree, n);
            long cc = (long)(tree.leafCount - 2) * (tree.leafCount - 3) / 2;
            if (cc <= 0) continue;
            for (int a = 0; a < n; a++) {
                int fa = td.firstOcc[a];
                if (fa < 0) continue;
                for (int b = a + 1; b < n; b++) {
                    int fb = td.firstOcc[b];
                    if (fb < 0) continue;
                    int lo = Math.min(fa, fb), hi = Math.max(fa, fb);
                    int width = hi - lo + 1;
                    int lvl = 31 - Integer.numberOfLeadingZeros(width);
                    int lo2 = hi - (1 << lvl) + 1;
                    int posL = td.sparseArgmin[lvl][lo];
                    int posR = td.sparseArgmin[lvl][lo2];
                    int pos = (td.depths[posL] <= td.depths[posR]) ? posL : posR;

                    int aS, bS;
                    double aF, bF;
                    if (fa <= fb) {
                        aS = td.eulerLeftChildS[pos];   aF = td.eulerLeftChildF[pos];
                        bS = td.eulerRightChildS[pos];  bF = td.eulerRightChildF[pos];
                    } else {
                        aS = td.eulerRightChildS[pos];  aF = td.eulerRightChildF[pos];
                        bS = td.eulerLeftChildS[pos];   bF = td.eulerLeftChildF[pos];
                    }
                    long z = tree.leafCount - aS - bS;
                    double twoQD = (td.eulerF[fa] - aF) + (td.eulerF[fb] - bF)
                        + (double)((long)(aS - 1) * z) + (double)((long)(bS - 1) * z);
                    double sameSide = (double)cc - twoQD * 0.5;
                    num[a * n + b] += sameSide; num[b * n + a] += sameSide;
                    den[a * n + b] += cc;       den[b * n + a] += cc;
                }
            }
        }

        SimilarityMatrix reference = SimilarityMatrixBuilder.buildCPU(trees, n);
        long checks = 0;
        for (int a = 0; a < n; a++) {
            for (int b = a + 1; b < n; b++) {
                double actual = den[a * n + b] > 0.0 ? num[a * n + b] / den[a * n + b] : 0.0;
                double expected = reference.getSim(a, b);
                if (Math.abs(actual - expected) > 1e-9) {
                    throw new AssertionError("compact similarity mismatch: pair=(" + a + "," + b
                        + ") expected=" + expected + " actual=" + actual);
                }
                checks++;
            }
        }
        return checks;
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            throw new IllegalArgumentException("Pass one or more Newick gene-tree files");
        }
        long treesChecked = 0;
        long intervalsChecked = 0;
        long matrixCellsChecked = 0;
        for (String path : args) {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> trees = TreeParser.parseGeneTrees(path, registry);
            for (Tree tree : trees) {
                try {
                    intervalsChecked += verifyTree(tree, registry.size());
                } catch (RuntimeException | AssertionError e) {
                    throw new AssertionError("Failure in " + path + " treeIndex="
                        + tree.treeIndex + " leaves=" + tree.leafCount, e);
                }
                treesChecked++;
            }
            matrixCellsChecked += verifyCompactSimilarity(trees, registry.size());
        }
        System.out.printf("Similarity argmin RMQ: PASS (%d trees, %d exact interval checks, "
                + "%d compact-vs-CPU matrix cells)%n",
            treesChecked, intervalsChecked, matrixCellsChecked);
    }
}
