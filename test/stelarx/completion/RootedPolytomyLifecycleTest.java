package stelarx.completion;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.ArrayList;
import java.util.List;

/** Invariant tests spanning native polytomy parsing, refinement, and completion. */
public final class RootedPolytomyLifecycleTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass one rooted-polytomy fixture");

        Parsed keptSerial = parse(args[0], true, 1);
        Parsed refinedSerial = parse(args[0], false, 1);
        Parsed refinedParallel = parse(args[0], false,
            Math.min(4, Runtime.getRuntime().availableProcessors()));
        check(keptSerial.trees.size() == refinedSerial.trees.size(), "tree count changed");

        int nativePolytomies = 0;
        for (int i = 0; i < keptSerial.trees.size(); i++) {
            Tree kept = keptSerial.trees.get(i);
            Tree serial = refinedSerial.trees.get(i);
            Tree parallel = refinedParallel.trees.get(i);
            validateRoot(kept);
            validateRoot(serial);
            validateRoot(parallel);
            if (kept.hasPolytomy) nativePolytomies++;
            check(!serial.hasPolytomy, "serial refinement retained a polytomy");
            check(!parallel.hasPolytomy, "parallel refinement retained a polytomy");
            check(serial.toNewick(refinedSerial.registry).equals(
                    parallel.toNewick(refinedParallel.registry)),
                "serial/parallel deterministic refinement mismatch at tree " + i);
        }
        check(nativePolytomies > 0, "fixture did not exercise a native internal polytomy");

        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        try {
            int n = keptSerial.registry.size();
            System.clearProperty("stelarx.similarity.forcePacked");
            SimilarityMatrix denseMatrix = SimilarityMatrixBuilder.buildCPU(keptSerial.trees, n);
            List<Tree> dense = TreeCompleter.completeAll(keptSerial.trees, denseMatrix, n);

            System.setProperty("stelarx.similarity.forcePacked", "true");
            SimilarityMatrix packedMatrix = SimilarityMatrixBuilder.buildCPU(keptSerial.trees, n);
            check(packedMatrix.isPacked(), "forced packed matrix was not selected");
            List<Tree> packed = TreeCompleter.completeAll(keptSerial.trees, packedMatrix, n);

            long triplesChecked = 0;
            for (int i = 0; i < keptSerial.trees.size(); i++) {
                Tree before = keptSerial.trees.get(i);
                Tree after = dense.get(i);
                Tree packedAfter = packed.get(i);
                validateComplete(after, n);
                validateComplete(packedAfter, n);
                check(after.toNewick(keptSerial.registry).equals(
                        packedAfter.toNewick(keptSerial.registry)),
                    "dense/packed completion mismatch at tree " + i);
                for (int a = 0; a < n; a++) {
                    if (before.positionMap[a] < 0) continue;
                    for (int b = a + 1; b < n; b++) {
                        if (before.positionMap[b] < 0) continue;
                        for (int c = b + 1; c < n; c++) {
                            if (before.positionMap[c] < 0) continue;
                            int expected = resolvedPair(before, before.root, a, b, c, n);
                            int actual = resolvedPair(after, after.root, a, b, c, n);
                            check(expected == actual, "completion changed rooted triple "
                                + a + "," + b + "," + c + " in tree " + i);
                            triplesChecked++;
                        }
                    }
                }
            }
            System.out.println("Rooted polytomy lifecycle: PASS (" + keptSerial.trees.size()
                + " trees, " + triplesChecked + " preserved triples)");
        } finally {
            System.clearProperty("stelarx.similarity.forcePacked");
            Threading.shutdown();
        }
    }

    private static Parsed parse(String path, boolean keep, int threads) throws Exception {
        TaxonRegistry registry = new TaxonRegistry();
        Threading.start(threads);
        try {
            return new Parsed(TreeParser.parseGeneTrees(path, registry, keep), registry);
        } finally {
            Threading.shutdown();
        }
    }

    private static void validateRoot(Tree tree) {
        check(tree.root != null && !tree.root.isLeaf(), "missing internal root");
        check(!tree.root.isPolytomous(), "supplied root ceased to be binary");
        check(tree.root.left != null && tree.root.right != null, "root child missing");
    }

    private static void validateComplete(Tree tree, int n) {
        validateRoot(tree);
        check(tree.isComplete && tree.leafCount == n, "completion did not restore all taxa");
        boolean[] seen = new boolean[n];
        validateNode(tree.root, null, seen);
        for (int taxon = 0; taxon < n; taxon++) {
            check(seen[taxon], "completed tree missing taxon " + taxon);
            check(tree.positionMap[taxon] >= 0, "completed position map missing taxon " + taxon);
        }
    }

    private static void validateNode(TreeNode node, TreeNode parent, boolean[] seen) {
        check(node.parent == parent, "broken parent pointer");
        if (node.isLeaf()) {
            check(node.taxonId >= 0 && node.taxonId < seen.length, "invalid taxon ID");
            check(!seen[node.taxonId], "duplicate taxon ID");
            seen[node.taxonId] = true;
            return;
        }
        for (TreeNode child : children(node)) validateNode(child, node, seen);
    }

    private static int resolvedPair(Tree tree, TreeNode node,
                                    int a, int b, int c, int n) {
        List<TreeNode> children = children(node);
        int ia = childContaining(tree, children, a);
        int ib = childContaining(tree, children, b);
        int ic = childContaining(tree, children, c);
        if (ia == ib && ia != ic) return a * n + b;
        if (ia == ic && ia != ib) return a * n + c;
        if (ib == ic && ib != ia) return b * n + c;
        if (ia == ib && ib == ic) return resolvedPair(tree, children.get(ia), a, b, c, n);
        return -1;
    }

    private static int childContaining(Tree tree, List<TreeNode> children, int taxon) {
        int position = tree.positionMap[taxon];
        for (int i = 0; i < children.size(); i++) {
            TreeNode child = children.get(i);
            if (position >= child.rangeStart && position < child.rangeEnd) return i;
        }
        throw new AssertionError("taxon not found below internal node");
    }

    private static List<TreeNode> children(TreeNode node) {
        List<TreeNode> result = new ArrayList<>();
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) result.add(child);
        } else {
            result.add(node.left);
            result.add(node.right);
        }
        return result;
    }

    private record Parsed(List<Tree> trees, TaxonRegistry registry) {}

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
