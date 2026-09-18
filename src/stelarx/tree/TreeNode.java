package stelarx.tree;

/**
 * Node in a rooted gene tree (binary or polytomous).
 * After parsing, every node has [rangeStart, rangeEnd) covering its leaf positions
 * in the tree's left-to-right postorder array.
 *
 * Polytomy representation:
 *   - Binary internal node: {@code left}/{@code right} set, {@code children == null}.
 *   - Polytomous internal node (k ≥ 3 children): {@code children} is a length-k array
 *     in left-to-right order; {@code left == children[0]} and {@code right == children[k-1]}
 *     so that {@code rangeStart}/{@code rangeEnd} and {@code isLeaf()} keep working.
 * The {@code children} array is the canonical child list when {@code isPolytomous()};
 * for binary nodes iterate {@code left}/{@code right} as before.
 */
public class TreeNode {
    public TreeNode left, right, parent;

    /** Non-null (length ≥ 3) for polytomous internal nodes; null otherwise. */
    public TreeNode[] children;

    /** Taxon ID for leaves; -1 for internal nodes. */
    public int taxonId = -1;

    /**
     * Half-open range [rangeStart, rangeEnd) in the tree's postorder array.
     * For a leaf: rangeEnd = rangeStart + 1.
     * For an internal node: spans entire descendant leaf range.
     */
    public int rangeStart = -1;
    public int rangeEnd   = -1;

    public boolean isLeaf() { return left == null; }
    public boolean isRoot() { return parent == null; }
    /** True for an internal node with k ≥ 3 children (its {@code children} array is set). */
    public boolean isPolytomous() { return children != null; }
    public int rangeSize()  { return rangeEnd - rangeStart; }

    /**
     * The other child of our parent (null if we are root).
     * NOTE: only well-defined when the parent is BINARY.  Must not be called for a
     * child of a polytomous parent (no unique sibling) — callers guard with
     * {@code !parent.isPolytomous()} (see DPTable §3.7).
     */
    public TreeNode getSibling() {
        if (parent == null) return null;
        return (parent.left == this) ? parent.right : parent.left;
    }
}
