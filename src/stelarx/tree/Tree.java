package stelarx.tree;

import stelarx.taxon.TaxonRegistry;

/**
 * A parsed rooted binary gene tree with postorder (L-to-R leaf) array
 * and inverse position map.
 */
public class Tree {
    /** Index in the gene tree list (0..k-1). */
    public final int treeIndex;

    /** Root of the tree. */
    public final TreeNode root;

    /**
     * postorderArray[pos] = taxon ID at left-to-right position pos.
     * Length = leafCount.
     */
    public final int[] postorderArray;

    /**
     * positionMap[taxonId] = position in postorderArray (-1 if absent).
     * Length = total taxa count n.
     */
    public final int[] positionMap;

    /** Number of leaves in this tree. */
    public final int leafCount;

    /** True when this tree contains all n taxa. */
    public final boolean isComplete;

    /** True when at least one internal node has three or more rooted children. */
    public final boolean hasPolytomy;

    public Tree(int treeIndex, TreeNode root,
                int[] postorderArray, int[] positionMap,
                int leafCount, int totalTaxa) {
        this(treeIndex, root, postorderArray, positionMap, leafCount, totalTaxa, false);
    }

    public Tree(int treeIndex, TreeNode root,
                int[] postorderArray, int[] positionMap,
                int leafCount, int totalTaxa, boolean hasPolytomy) {
        this.treeIndex = treeIndex;
        this.root = root;
        this.postorderArray = postorderArray;
        this.positionMap = positionMap;
        this.leafCount = leafCount;
        this.isComplete = (leafCount == totalTaxa);
        this.hasPolytomy = hasPolytomy;
    }

    /** Reconstruct Newick string (no branch lengths). */
    public String toNewick(TaxonRegistry reg) {
        return (hasPolytomy ? nodeToNewickPolytomy(root, reg) : nodeToNewick(root, reg)) + ";";
    }

    private String nodeToNewick(TreeNode n, TaxonRegistry reg) {
        if (n.isLeaf()) return NewickLabel.encode(reg.getName(n.taxonId));
        return "(" + nodeToNewick(n.left, reg) + "," + nodeToNewick(n.right, reg) + ")";
    }

    private String nodeToNewickPolytomy(TreeNode n, TaxonRegistry reg) {
        if (n.isLeaf()) return NewickLabel.encode(reg.getName(n.taxonId));
        if (!n.isPolytomous()) {
            return "(" + nodeToNewickPolytomy(n.left, reg) + ","
                + nodeToNewickPolytomy(n.right, reg) + ")";
        }
        StringBuilder out = new StringBuilder("(");
        for (int i = 0; i < n.children.length; i++) {
            if (i > 0) out.append(',');
            out.append(nodeToNewickPolytomy(n.children[i], reg));
        }
        return out.append(')').toString();
    }
}
