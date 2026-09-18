package stelarx.completion;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.List;

/** Regression coverage for completion of incomplete, natively polytomous trees. */
public final class TreeCompleterPolytomyTest {
    private static final String[] EXPECTED = {
        "(F,((E,D),(A,G,(B,C))));",
        "(F,(G,(E,(A,(C,D,B)))));",
        "(G,(F,(D,(B,(C,E,A)))));",
        "(A,(C,(D,(F,(E,G,B)))));"
    };

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass the incomplete-polytomy fixture");
        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        try {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> trees = TreeParser.parseGeneTrees(args[0], registry, true);
            for (Tree tree : trees) {
                check(tree.hasPolytomy, "input tree " + tree.treeIndex + " lost its polytomy marker");
                check(!tree.isComplete, "fixture tree " + tree.treeIndex + " must be incomplete");
            }

            int n = registry.size();
            System.clearProperty("stelarx.similarity.forcePacked");
            SimilarityMatrix denseMatrix = SimilarityMatrixBuilder.buildCPU(trees, n);
            List<Tree> dense = TreeCompleter.completeAll(trees, denseMatrix, n);

            System.setProperty("stelarx.similarity.forcePacked", "true");
            SimilarityMatrix packedMatrix;
            try {
                packedMatrix = SimilarityMatrixBuilder.buildCPU(trees, n);
            } finally {
                System.clearProperty("stelarx.similarity.forcePacked");
            }
            check(packedMatrix.isPacked(), "packed completion path was not selected");
            List<Tree> packed = TreeCompleter.completeAll(trees, packedMatrix, n);

            for (int i = 0; i < trees.size(); i++) {
                validateCompleted(dense.get(i), n);
                validateCompleted(packed.get(i), n);
                String denseNewick = dense.get(i).toNewick(registry);
                String packedNewick = packed.get(i).toNewick(registry);
                check(denseNewick.equals(packedNewick),
                    "dense/packed completion mismatch for tree " + i);
                check(EXPECTED[i].equals(denseNewick),
                    "completion topology changed for tree " + i + ": " + denseNewick);
            }
            System.out.println("Incomplete polytomy completion: PASS");
        } finally {
            Threading.shutdown();
            System.clearProperty("stelarx.similarity.forcePacked");
        }
    }

    private static void validateCompleted(Tree tree, int n) {
        check(tree.isComplete, "completed tree " + tree.treeIndex + " is still incomplete");
        check(tree.leafCount == n, "completed tree " + tree.treeIndex + " has wrong leaf count");
        check(tree.hasPolytomy, "completed tree " + tree.treeIndex + " lost all polytomies");
        boolean[] seen = new boolean[n];
        int[] position = {0};
        int[] polytomies = {0};
        validateNode(tree.root, null, seen, position, polytomies);
        check(position[0] == n, "tree traversal did not reach every taxon");
        check(polytomies[0] > 0, "tree contains no canonical polytomy node");
        for (int taxon = 0; taxon < n; taxon++) {
            check(seen[taxon], "missing taxon " + taxon);
            check(tree.positionMap[taxon] >= 0, "position map missing taxon " + taxon);
        }
    }

    private static void validateNode(TreeNode node, TreeNode parent, boolean[] seen,
                                     int[] position, int[] polytomies) {
        check(node.parent == parent, "broken parent pointer");
        if (node.isLeaf()) {
            check(node.taxonId >= 0 && node.taxonId < seen.length, "invalid leaf taxon ID");
            check(!seen[node.taxonId], "duplicate leaf taxon " + node.taxonId);
            seen[node.taxonId] = true;
            check(node.rangeStart == position[0] && node.rangeEnd == position[0] + 1,
                "invalid leaf range");
            position[0]++;
            return;
        }

        int start = position[0];
        if (node.isPolytomous()) {
            polytomies[0]++;
            check(node.children.length >= 3, "polytomy has fewer than three children");
            check(node.left == node.children[0], "left alias is not first polytomy child");
            check(node.right == node.children[node.children.length - 1],
                "right alias is not last polytomy child");
            for (TreeNode child : node.children) validateNode(child, node, seen, position, polytomies);
        } else {
            check(node.left != null && node.right != null, "binary node is missing a child");
            validateNode(node.left, node, seen, position, polytomies);
            validateNode(node.right, node, seen, position, polytomies);
        }
        check(node.rangeStart == start && node.rangeEnd == position[0], "invalid internal range");
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
