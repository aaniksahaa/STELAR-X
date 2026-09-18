package stelarx.tree;

import stelarx.taxon.TaxonRegistry;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Topology-preserving induced-subtree restriction with unary-node suppression. */
public final class TreeRestrictor {
    private TreeRestrictor() {}

    /**
     * Restrict gene trees in place at the node level and rebuild their compact
     * arrays against {@code targetRegistry}. Trees with fewer than three retained
     * taxa are dropped because they contain no rooted triplet and therefore contribute
     * exactly zero to every score.
     */
    public static GeneResult restrictGeneTrees(List<Tree> sourceTrees,
                                               TaxonRegistry sourceRegistry,
                                               TaxonRegistry targetRegistry) {
        Map<String, Integer> targetIds = targetIdMap(targetRegistry);
        List<Tree> restricted = new ArrayList<>(sourceTrees.size());
        int dropped = 0;
        for (Tree source : sourceTrees) {
            Tree tree = restrictOne(source, sourceRegistry, targetRegistry,
                targetIds, restricted.size());
            if (tree == null || tree.leafCount < 3) {
                dropped++;
            } else {
                restricted.add(tree);
            }
        }
        if (restricted.isEmpty()) {
            throw new IllegalArgumentException(
                "No gene tree retains at least three taxa after applying the taxa file");
        }
        return new GeneResult(restricted, dropped);
    }

    /** Restrict the supplied species tree and require every target taxon exactly once. */
    public static Tree restrictSpeciesTree(Tree source,
                                           TaxonRegistry sourceRegistry,
                                           TaxonRegistry targetRegistry) {
        Map<String, Integer> targetIds = targetIdMap(targetRegistry);
        Tree restricted = restrictOne(source, sourceRegistry, targetRegistry, targetIds, 0);
        if (restricted == null || restricted.leafCount < 3) {
            throw new IllegalArgumentException(
                "Species tree retains fewer than three scoring taxa after applying the taxa file");
        }
        if (restricted.leafCount != targetRegistry.size()) {
            throw new IllegalArgumentException("Restricted species tree contains "
                + restricted.leafCount + " of " + targetRegistry.size()
                + " effective scoring taxa");
        }
        return restricted;
    }

    private static Tree restrictOne(Tree source,
                                    TaxonRegistry sourceRegistry,
                                    TaxonRegistry targetRegistry,
                                    Map<String, Integer> targetIds,
                                    int newTreeIndex) {
        TreeNode root = prune(source.root, sourceRegistry, targetIds);
        if (root == null) return null;
        root.parent = null;

        int n = targetRegistry.size();
        int[] order = new int[n];
        int[] counter = {0};
        assignRanges(root, order, counter);
        int leafCount = counter[0];
        order = Arrays.copyOf(order, leafCount);

        int[] positions = new int[n];
        Arrays.fill(positions, -1);
        for (int pos = 0; pos < leafCount; pos++) {
            int taxon = order[pos];
            if (positions[taxon] >= 0) {
                throw new IllegalArgumentException("Tree " + source.treeIndex
                    + " contains duplicate taxon after restriction: "
                    + targetRegistry.getName(taxon));
            }
            positions[taxon] = pos;
        }
        return new Tree(newTreeIndex, root, order, positions, leafCount, n,
            containsPolytomy(root));
    }

    /** Returns null for an entirely removed clade and suppresses unary nodes. */
    private static TreeNode prune(TreeNode node,
                                  TaxonRegistry sourceRegistry,
                                  Map<String, Integer> targetIds) {
        if (node.isLeaf()) {
            Integer mapped = targetIds.get(sourceRegistry.getName(node.taxonId));
            if (mapped == null) return null;
            node.taxonId = mapped;
            node.parent = null;
            return node;
        }

        List<TreeNode> kept = new ArrayList<>(node.isPolytomous()
            ? node.children.length : 2);
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) {
                TreeNode retained = prune(child, sourceRegistry, targetIds);
                if (retained != null) kept.add(retained);
            }
        } else {
            TreeNode left = prune(node.left, sourceRegistry, targetIds);
            TreeNode right = prune(node.right, sourceRegistry, targetIds);
            if (left != null) kept.add(left);
            if (right != null) kept.add(right);
        }

        if (kept.isEmpty()) return null;
        if (kept.size() == 1) {
            TreeNode only = kept.get(0);
            only.parent = null;
            return only;
        }

        node.taxonId = -1;
        if (kept.size() == 2) {
            node.children = null;
            node.left = kept.get(0);
            node.right = kept.get(1);
        } else {
            node.children = kept.toArray(TreeNode[]::new);
            node.left = node.children[0];
            node.right = node.children[node.children.length - 1];
        }
        for (TreeNode child : kept) child.parent = node;
        return node;
    }

    private static void assignRanges(TreeNode node, int[] order, int[] counter) {
        if (node.isLeaf()) {
            node.rangeStart = counter[0];
            node.rangeEnd = counter[0] + 1;
            order[counter[0]++] = node.taxonId;
            return;
        }
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) assignRanges(child, order, counter);
        } else {
            assignRanges(node.left, order, counter);
            assignRanges(node.right, order, counter);
        }
        node.rangeStart = node.left.rangeStart;
        node.rangeEnd = node.right.rangeEnd;
    }

    private static boolean containsPolytomy(TreeNode node) {
        if (node.isLeaf()) return false;
        if (node.isPolytomous()) return true;
        return containsPolytomy(node.left) || containsPolytomy(node.right);
    }

    private static Map<String, Integer> targetIdMap(TaxonRegistry registry) {
        if (!registry.isLocked()) {
            throw new IllegalArgumentException("Target taxon registry must be locked");
        }
        Map<String, Integer> ids = new HashMap<>(Math.max(16, registry.size() * 2));
        for (int i = 0; i < registry.size(); i++) ids.put(registry.getName(i), i);
        return ids;
    }

    public record GeneResult(List<Tree> trees, int droppedTreeCount) {}
}
