package stelarx.completion;

import stelarx.Logging;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.ProgressBar;
import stelarx.util.Threading;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Inserts missing taxa into incomplete rooted gene trees while preserving the
 * supplied root and every relationship among the taxa already present.
 *
 * Algorithm for inserting taxon x into tree T:
 *   1. Find anchor a = closest taxon currently present in T (via sortedRows).
 *   2. Replace leaf a by a new binary parent whose children are a and x.
 *   3. Update the present-taxon and leaf-node maps.
 *
 * This never reverses an edge. Deleting inserted leaves and suppressing their
 * degree-one parents therefore recovers exactly the supplied rooted tree.
 *
 * Missing taxa are processed in ascending ID order for each tree.
 * Different trees are completed in parallel using the thread pool.
 *
 * The similarity ranking and packed/dense matrix implementations are retained
 * from STELAR-X; only the topology mutation is root-aware.
 */
public class TreeCompleter {

    /**
     * Complete all incomplete trees in the list using the four-point algorithm.
     *
     * @param trees      gene trees (may mix complete and incomplete)
     * @param sim        flat n×n similarity matrix (sim[a*n+b] ∈ [0,1]; used for four-point)
     * @param dist       flat n×n distance matrix (dist = 1-sim; used for sortedRows)
     * @param n          total taxon count
     * @return new list where every tree is complete; already-complete trees pass through
     */
    public static List<Tree> completeAll(List<Tree> trees, double[] sim, double[] dist, int n) {
        List<Integer> incomplete = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            Tree tree = trees.get(i);
            if (!tree.isComplete) {
                incomplete.add(i);
            }
        }

        if (incomplete.isEmpty()) return trees;

        Logging.info("Tree completion: %d/%d trees incomplete", incomplete.size(), trees.size());

        // Build sortedRows once — shared read-only across all tree completions.
        // sortedRows[x*n + rank] = taxon ID of x's rank-th nearest neighbor (ascending dist).
        int[] sortedRows = SortedRowsBuilder.buildCPU(dist, n);

        // Mutable result array; complete trees pass through unchanged.
        Tree[] result = trees.toArray(new Tree[0]);

        ProgressBar bar   = new ProgressBar("Completing incomplete gene trees", incomplete.size());
        AtomicInteger cnt = new AtomicInteger(0);

        PolyMatrix matrix = new DensePolyMatrix(sim, sortedRows, n);
        Threading.processParallel(incomplete, idx -> {
            result[idx] = completeTreeRootPreserving(trees.get(idx), matrix, n);
            bar.update(cnt.incrementAndGet());
        });
        bar.done();

        return Arrays.asList(result);
    }

    /** Exact large-matrix entry point; the original dense path remains unchanged. */
    public static List<Tree> completeAll(List<Tree> trees, SimilarityMatrix sim, int n) {
        if (!sim.isPacked()) return completeAll(trees, sim.sim, sim.dist, n);

        List<Integer> incomplete = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            Tree tree = trees.get(i);
            if (!tree.isComplete) {
                incomplete.add(i);
            }
        }
        if (incomplete.isEmpty()) return trees;

        Logging.info("Tree completion: %d/%d trees incomplete (segmented exact matrix)",
            incomplete.size(), trees.size());
        int[][] sortedRows = SortedRowsBuilder.buildPackedCPU(sim);
        Tree[] result = trees.toArray(new Tree[0]);
        ProgressBar bar = new ProgressBar("Completing incomplete gene trees", incomplete.size());
        AtomicInteger cnt = new AtomicInteger(0);
        PolyMatrix matrix = new PackedPolyMatrix(sim, sortedRows);
        Threading.processParallel(incomplete, idx -> {
            result[idx] = completeTreeRootPreserving(trees.get(idx), matrix, n);
            bar.update(cnt.incrementAndGet());
        });
        bar.done();
        return Arrays.asList(result);
    }

    /**
     * Root-aware completion used by STELAR-X. A missing taxon is attached as the
     * sibling of its most similar present taxon. Replacing that leaf by a new
     * binary parent never reverses an edge and therefore preserves the input root.
     */
    private static Tree completeTreeRootPreserving(Tree tree, PolyMatrix matrix, int n) {
        boolean[] inTree = new boolean[n];
        TreeNode[] taxonNode = new TreeNode[n];
        TreeNode root = tree.hasPolytomy
            ? deepCopyNodesPolytomy(tree.root, null, taxonNode)
            : deepCopyNodes(tree.root, null, taxonNode);
        for (int i = 0; i < n; i++) if (tree.positionMap[i] != -1) inTree[i] = true;

        for (int x = 0; x < n; x++) {
            if (inTree[x]) continue;
            int anchor = matrix.findAnchor(x, inTree);
            TreeNode anchorLeaf = taxonNode[anchor];
            if (anchorLeaf == null) {
                throw new IllegalStateException("Missing completion anchor node for taxon " + anchor);
            }
            TreeNode oldParent = anchorLeaf.parent;
            TreeNode newLeaf = new TreeNode();
            newLeaf.taxonId = x;
            TreeNode joined = new TreeNode();
            setBinaryChildren(joined, anchorLeaf, newLeaf);
            anchorLeaf.parent = joined;
            newLeaf.parent = joined;
            joined.parent = oldParent;
            if (oldParent == null) root = joined;
            else replaceChild(oldParent, anchorLeaf, joined);
            inTree[x] = true;
            taxonNode[x] = newLeaf;
        }

        return tree.hasPolytomy
            ? rebuildTreePolytomy(tree.treeIndex, root, n)
            : rebuildTree(tree.treeIndex, root, n);
    }

    private static Tree completeTreeFourPointPacked(Tree tree, SimilarityMatrix sim,
                                                     int[][] sortedRows, int n) {
        boolean[] inTree = new boolean[n];
        TreeNode[] taxonNode = new TreeNode[n];
        TreeNode root = deepCopyNodes(tree.root, null, taxonNode);
        for (int i = 0; i < n; i++) if (tree.positionMap[i] != -1) inTree[i] = true;
        root = preprocessReroot(root, tree.leafCount);

        for (int x = 0; x < n; x++) {
            if (inTree[x]) continue;
            int anchor = findAnchorPacked(x, inTree, sortedRows);
            root = rerootAtLeafEdge(taxonNode[anchor], root);
            TreeNode start = root.right;
            int c1rep = -1, c2rep = -1;
            TreeNode c1 = null, c2 = null;
            while (!start.isLeaf()) {
                c1 = start.left;
                c2 = start.right;
                if (c1rep == -1) c1rep = leftmostTaxon(c1);
                if (c2rep == -1) c2rep = leftmostTaxon(c2);
                int better = fourPointBetterSidePacked(x, anchor, c1rep, c2rep, sim);
                if (better == anchor) break;
                if (better == c1rep) {
                    start = c1;
                    c2rep = -1;
                } else {
                    start = c2;
                    c1rep = c2rep;
                    c2rep = -1;
                }
            }
            TreeNode newLeaf = insertTaxon(x, start, c1, c2);
            inTree[x] = true;
            taxonNode[x] = newLeaf;
        }
        return rebuildTree(tree.treeIndex, root, n);
    }

    private static int findAnchorPacked(int x, boolean[] inTree, int[][] sortedRows) {
        int[] row = sortedRows[x];
        for (int candidate : row) {
            if (candidate != x && inTree[candidate]) return candidate;
        }
        throw new RuntimeException("No anchor found for taxon " + x + " — inTree array may be empty");
    }

    private static int fourPointBetterSidePacked(int x, int a, int b, int c,
                                                  SimilarityMatrix sim) {
        double xa = sim.getSim(x, a);
        double xb = sim.getSim(x, b);
        double xc = sim.getSim(x, c);
        double ab = sim.getSim(a, b);
        double ac = sim.getSim(a, c);
        double bc = sim.getSim(b, c);
        double ascore = xa + bc - (xb + ac);
        double bscore = xb + ac - (xa + bc);
        double cscore = xc + ab - (xb + ac);
        return ascore >= bscore
            ? (ascore >= cscore ? a : c)
            : (bscore >= cscore ? b : c);
    }

    // ── Polytomous-tree completion ──────────────────────────────────────────

    /**
     * Completion path for trees containing native polytomies.  Binary trees do
     * not enter this path, preserving their established allocation and traversal
     * behaviour.  Child-list mutations below mirror STITree's remove/adopt/append
     * semantics used by ASTRAL-MP's completion routine without discarding any
     * unresolved arms.
     */
    private static Tree completeTreeFourPointPolytomy(Tree tree, PolyMatrix matrix, int n) {
        boolean[] inTree = new boolean[n];
        TreeNode[] taxonNode = new TreeNode[n];
        TreeNode root = deepCopyNodesPolytomy(tree.root, null, taxonNode);
        for (int i = 0; i < n; i++) if (tree.positionMap[i] != -1) inTree[i] = true;
        root = preprocessRerootPolytomy(root, tree.leafCount);

        for (int x = 0; x < n; x++) {
            if (inTree[x]) continue;
            int anchor = matrix.findAnchor(x, inTree);
            TreeNode anchorLeaf = taxonNode[anchor];
            if (anchorLeaf == null) {
                throw new IllegalStateException("Tree " + tree.treeIndex + ": anchor taxon "
                    + anchor + " is marked present but has no leaf node");
            }
            root = rerootAtEdgePolytomy(anchorLeaf, root);
            TreeNode start = root.right;
            int c1rep = -1, c2rep = -1;
            TreeNode c1 = null, c2 = null;

            while (!start.isLeaf()) {
                c1 = childAt(start, 0);
                c2 = childAt(start, 1);
                if (c1rep == -1) c1rep = leftmostTaxon(c1);
                if (c2rep == -1) c2rep = leftmostTaxon(c2);
                int better = matrix.betterSide(x, anchor, c1rep, c2rep);
                if (better == anchor) break;
                if (better == c1rep) {
                    start = c1;
                    c2rep = -1;
                } else {
                    start = c2;
                    c1rep = c2rep;
                    c2rep = -1;
                }
            }

            TreeNode newLeaf = insertTaxonPolytomy(x, start, c1, c2);
            inTree[x] = true;
            taxonNode[x] = newLeaf;
        }
        return rebuildTreePolytomy(tree.treeIndex, root, n);
    }

    private interface PolyMatrix {
        int findAnchor(int x, boolean[] inTree);
        int betterSide(int x, int a, int b, int c);
    }

    private static final class DensePolyMatrix implements PolyMatrix {
        private final double[] sim;
        private final int[] sortedRows;
        private final int n;

        DensePolyMatrix(double[] sim, int[] sortedRows, int n) {
            this.sim = sim;
            this.sortedRows = sortedRows;
            this.n = n;
        }

        public int findAnchor(int x, boolean[] inTree) {
            return TreeCompleter.findAnchor(x, inTree, sortedRows, n);
        }

        public int betterSide(int x, int a, int b, int c) {
            return fourPointBetterSide(x, a, b, c, sim, n);
        }
    }

    private static final class PackedPolyMatrix implements PolyMatrix {
        private final SimilarityMatrix sim;
        private final int[][] sortedRows;

        PackedPolyMatrix(SimilarityMatrix sim, int[][] sortedRows) {
            this.sim = sim;
            this.sortedRows = sortedRows;
        }

        public int findAnchor(int x, boolean[] inTree) {
            return findAnchorPacked(x, inTree, sortedRows);
        }

        public int betterSide(int x, int a, int b, int c) {
            return fourPointBetterSidePacked(x, a, b, c, sim);
        }
    }

    // ── Per-tree completion ───────────────────────────────────────────────────

    /**
     * Insert every taxon missing from this tree using the four-point algorithm
     * and return a rebuilt Tree.
     *
     * Mirrors WQDataCollection.getCompleteTree() lines 261–341.
     */
    private static Tree completeTreeFourPoint(Tree tree, double[] sim,
                                              int[] sortedRows, int n) {
        // --- Setup ---
        boolean[]  inTree    = new boolean[n];
        TreeNode[] taxonNode = new TreeNode[n];

        // Deep-copy tree nodes before any mutation; also populates taxonNode.
        TreeNode root = deepCopyNodes(tree.root, null, taxonNode);

        // Initialise inTree from the original positionMap.
        for (int i = 0; i < n; i++) {
            if (tree.positionMap[i] != -1) inTree[i] = true;
        }

        // --- Preprocessing reroot (mirrors WQDataCollection.reroot()) ---
        // ASTRAL-MP rereroots each gene tree at the most balanced internal node
        // (leaf count closest to leafCount/2) before the insertion loop begins,
        // and also moves direct leaf-children to the right slot so that
        // leftmostTaxon() always descends into an internal child first.
        // Both effects are deterministic and must be replicated to get identical
        // completed trees.
        root = preprocessReroot(root, tree.leafCount);

        // --- Insert each missing taxon in ascending ID order ---
        // (mirrors WQDataCollection loop: nextClearBit ascending)
        for (int x = 0; x < n; x++) {
            if (inTree[x]) continue;   // already present

            // Phase B: find anchor (closest in-tree taxon)
            // Mirrors AbstractMatrix.getClosestPresentTaxonId() lines 52–67.
            int anchor = findAnchor(x, inTree, sortedRows, n);

            // Phase C: physically reroot at edge (anchorLeaf, anchorLeaf.parent)
            // After rerooting: newRoot.left = anchorLeaf, newRoot.right = start.
            // Mirrors trc.rerootTreeAtNode(closestNode) + removeBinaryNodes in WQDataCollection.java:280.
            TreeNode anchorLeaf = taxonNode[anchor];
            root = rerootAtLeafEdge(anchorLeaf, root);

            // Phase D: four-point navigation
            // Mirrors the while(true) loop in WQDataCollection.java lines 290–323.
            TreeNode start = root.right;   // non-anchor child = "rest of tree"
            int c1rep = -1, c2rep = -1;
            TreeNode c1 = null, c2 = null;

            while (!start.isLeaf()) {
                // HINT for n-ary: c1 = children.get(0), c2 = children.get(1)
                c1 = start.left;
                c2 = start.right;

                if (c1rep == -1) c1rep = leftmostTaxon(c1);
                if (c2rep == -1) c2rep = leftmostTaxon(c2);

                int better = fourPointBetterSide(x, anchor, c1rep, c2rep, sim, n);

                if (better == anchor) {
                    // x groups with anchor → insert at this internal node.
                    break;
                } else if (better == c1rep) {
                    // Descend into left child.
                    // c1rep is still valid for the new left child (leftmost is preserved).
                    start = c1;
                    c2rep = -1;   // right side changes at next level
                } else {
                    // better == c2rep: descend into right child.
                    // c2's leftmost becomes the new c1rep.
                    start = c2;
                    c1rep = c2rep;
                    c2rep = -1;
                }
            }

            // Phase E: insert taxon x
            // c1/c2 hold the last-seen children (used for internal-node case).
            // Mirrors WQDataCollection.java lines 325–337.
            TreeNode newLeaf = insertTaxon(x, start, c1, c2);

            // Phase F: update membership so future iterations can use x as anchor.
            inTree[x]    = true;
            taxonNode[x] = newLeaf;
        }

        return rebuildTree(tree.treeIndex, root, n);
    }

    // ── Anchor finding ────────────────────────────────────────────────────────

    /**
     * Find the closest taxon to x that is currently in the tree.
     *
     * Scans sortedRows[x] in ascending distance order and returns the first
     * in-tree taxon found.
     *
     * Mirrors AbstractMatrix.getClosestPresentTaxonId() lines 52–67.
     * The original checks (missingId > other || presentBS.get(other)); here we
     * unify both conditions into inTree[candidate] which is true for both
     * originally-present and already-inserted taxa.
     */
    private static int findAnchor(int x, boolean[] inTree, int[] sortedRows, int n) {
        int base = x * n;
        for (int rank = 0; rank < n; rank++) {
            int candidate = sortedRows[base + rank];
            if (candidate != x && inTree[candidate]) {
                return candidate;
            }
        }
        throw new RuntimeException("No anchor found for taxon " + x
                + " — inTree array may be empty");
    }

    // ── Four-point score ──────────────────────────────────────────────────────

    /**
     * Determine which of {a (anchor), b (c1rep), c (c2rep)} taxon x groups with,
     * using the four-point condition on the similarity matrix.
     *
     * Scores (higher = better grouping):
     *   ascore = sim[x][a] + sim[b][c] − sim[x][b] − sim[a][c]   (x with anchor)
     *   bscore = sim[x][b] + sim[a][c] − sim[x][a] − sim[b][c]   (x with c1rep)
     *   cscore = sim[x][c] + sim[a][b] − sim[x][b] − sim[a][c]   (x with c2rep)
     *
     * Returns the taxon ID (a, b, or c) of the winning side.
     *
     * Mirrors SimilarityMatrix.getBetterSideByFourPoint() lines 45–58.
     */
    private static int fourPointBetterSide(int x, int a, int b, int c,
                                           double[] sim, int n) {
        double xa = sim[x * n + a];
        double xb = sim[x * n + b];
        double xc = sim[x * n + c];
        double ab = sim[a * n + b];
        double ac = sim[a * n + c];
        double bc = sim[b * n + c];

        double ascore = xa + bc - (xb + ac);
        double bscore = xb + ac - (xa + bc);
        double cscore = xc + ab - (xb + ac);

        // Mirrors the ternary in SimilarityMatrix.getBetterSideByFourPoint exactly.
        return ascore >= bscore
                ? (ascore >= cscore ? a : c)
                : (bscore >= cscore ? b : c);
    }

    // ── Physical rerooting ────────────────────────────────────────────────────

    /**
     * Physically reroot the tree at the edge between anchorLeaf and its parent.
     *
     * After this call:
     *   newRoot.left  = anchorLeaf
     *   newRoot.right = p1 (the former parent of anchorLeaf, now head of the
     *                   "start" subtree containing everything except anchorLeaf)
     *
     * Algorithm (path reversal):
     *   Collect path from anchorLeaf to old root: [anchor, p1, p2, ..., pk].
     *   Allocate newRoot with left=anchor, right=p1.
     *   Reverse parent→child edges along the path.
     *   Collapse the old root (it becomes a unary node after reversal and is
     *   spliced out so the tree stays binary).
     *
     * Mirrors WQDataCollection.java:280 (trc.rerootTreeAtNode + removeBinaryNodes).
     *
     * Complexity: O(depth) time and space.
     *
     * HINT for n-ary extension: in step 2, replace left/right slot with an entry
     * in node.children list. In the collapse step, remove the old-root entry from
     * path[k-1].children.
     *
     * @param anchorLeaf the leaf node that will become newRoot.left
     * @param oldRoot    the current root of the tree (needed only to detect path end)
     * @return the new root node
     */
    static TreeNode rerootAtLeafEdge(TreeNode anchorLeaf, TreeNode oldRoot) {
        // Build path from anchorLeaf up to (and including) old root.
        // path[0] = anchorLeaf, path[k] = oldRoot.
        List<TreeNode> path = new ArrayList<>();
        TreeNode cur = anchorLeaf;
        while (cur != null) {
            path.add(cur);
            if (cur == oldRoot) break;
            cur = cur.parent;
        }
        int k = path.size() - 1;   // index of old root

        // Allocate new root: left = anchor, right = p1.
        TreeNode newRoot = new TreeNode();
        TreeNode p1      = path.get(1);   // former parent of anchorLeaf
        newRoot.left     = anchorLeaf;
        newRoot.right    = p1;
        anchorLeaf.parent = newRoot;

        // Special case: anchorLeaf is a direct child of old root (path length = 2).
        // path = [anchor, root].  After newRoot.right = p1 = root, we need to
        // collapse the old root (now unary after anchor is stolen) into its
        // remaining child.
        if (k == 1) {
            // p1 == oldRoot.  Find the sibling of anchor under old root.
            TreeNode sib = (oldRoot.left == anchorLeaf) ? oldRoot.right : oldRoot.left;
            // Replace old root with sib directly.
            newRoot.right  = sib;
            sib.parent     = newRoot;
            // Old root is discarded.
            return newRoot;
        }

        // General case: path length >= 3 (anchor, p1, ..., pk=oldRoot).
        //
        // Step 2: reverse edges for nodes p1 .. pk-1 (indices 1 .. k-1) using
        // ASTRAL-MP's "remove + adopt to end" semantics:
        //   removeChild(childOnPath) then adoptChild(parentOnPath).
        // In LinkedList terms: the old parent ends up at the END of the
        // children list. For our binary representation that means:
        //   node.left  = otherChild  (the child not on the path)
        //   node.right = parentOnPath (former parent, now appended)
        // This is critical for replicating ASTRAL-MP because getLeftmostLeaf()
        // returns the post-order-first leaf, which depends on .left first.
        for (int i = 1; i <= k - 1; i++) {
            TreeNode node         = path.get(i);
            TreeNode childOnPath  = path.get(i - 1);   // toward anchor (keep as child)
            TreeNode parentOnPath = path.get(i + 1);   // old parent (becomes new child)

            TreeNode otherChild = (node.left == childOnPath) ? node.right : node.left;

            node.left  = otherChild;
            node.right = parentOnPath;

            // Fix parent pointers.
            node.parent         = (i == 1) ? newRoot : path.get(i - 1);
            parentOnPath.parent = node;
            // otherChild.parent is already `node`.
        }

        // Step 3: collapse old root pk using the same "remove + adopt to end"
        // semantics. After step 2:
        //   pathKm1.left  = (something not pk)
        //   pathKm1.right = pk
        // We want to remove pk from pathKm1 and adopt pk's remaining child at
        // the end. For binary: pathKm1.left stays, pathKm1.right = remainingChild.
        TreeNode pk             = path.get(k);   // old root
        TreeNode pathKm1        = path.get(k - 1);
        TreeNode remainingChild = (pk.left == pathKm1) ? pk.right : pk.left;

        // pathKm1's children after step 2 are [otherChild, pk]; replace pk → remainingChild.
        // Since pk is in the .right slot (per step 2), we set .right = remainingChild.
        pathKm1.right = remainingChild;
        remainingChild.parent = pathKm1;
        // pk is now unreachable and will be garbage-collected.

        return newRoot;
    }

    /** Polytomy-preserving counterpart of {@link #rerootAtLeafEdge}. */
    private static TreeNode rerootAtEdgePolytomy(TreeNode edgeNode, TreeNode oldRoot) {
        List<TreeNode> path = new ArrayList<>();
        TreeNode cur = edgeNode;
        while (cur != null) {
            path.add(cur);
            if (cur == oldRoot) break;
            cur = cur.parent;
        }
        if (path.size() < 2 || path.get(path.size() - 1) != oldRoot) {
            throw new IllegalStateException("Cannot reroot: target node is not below the current root");
        }

        int k = path.size() - 1;
        TreeNode newRoot = new TreeNode();
        TreeNode p1 = path.get(1);
        setBinaryChildren(newRoot, edgeNode, p1);
        edgeNode.parent = newRoot;

        if (k == 1) {
            requireBinaryRoot(oldRoot);
            TreeNode sibling = oldRoot.left == edgeNode ? oldRoot.right : oldRoot.left;
            setBinaryChildren(newRoot, edgeNode, sibling);
            edgeNode.parent = newRoot;
            sibling.parent = newRoot;
            return newRoot;
        }

        for (int i = 1; i <= k - 1; i++) {
            TreeNode node = path.get(i);
            TreeNode childOnPath = path.get(i - 1);
            TreeNode parentOnPath = path.get(i + 1);
            removeAndAppend(node, childOnPath, parentOnPath);
            node.parent = (i == 1) ? newRoot : path.get(i - 1);
            parentOnPath.parent = node;
        }

        TreeNode oldRootNode = path.get(k);
        TreeNode pathKm1 = path.get(k - 1);
        requireBinaryRoot(oldRootNode);
        TreeNode remaining = oldRootNode.left == pathKm1
            ? oldRootNode.right : oldRootNode.left;
        replaceChild(pathKm1, oldRootNode, remaining);
        remaining.parent = pathKm1;
        return newRoot;
    }

    // ── Insertion ─────────────────────────────────────────────────────────────

    /**
     * Insert taxon x at the position determined by the navigation loop.
     *
     * Two cases (mirrors WQDataCollection.java lines 325–337):
     *
     * Case 1 — stopped at a leaf (start.isLeaf()):
     *   Graft x as a sibling of start under a new internal node.
     *     start.parent → newInternal → {start, newLeaf(x)}
     *
     * Case 2 — stopped at an internal node (betterSide == anchor):
     *   x becomes a new direct child of start; c1 and c2 are wrapped under
     *   a new internal node as the other child.
     *     start → {newLeaf(x), newInternal → {c1, c2}}
     *
     * HINT for n-ary: in case 2, wrap ALL existing children under newInternal:
     *   newInternal.children = start.children; start.children = [newLeaf, newInternal].
     *
     * @param x     taxon ID to insert
     * @param start node where navigation stopped
     * @param c1    start.left as of the last loop iteration (used in case 2 only)
     * @param c2    start.right as of the last loop iteration (used in case 2 only)
     * @return the newly created leaf node for taxon x
     */
    private static TreeNode insertTaxon(int x, TreeNode start,
                                        TreeNode c1, TreeNode c2) {
        TreeNode newLeaf = new TreeNode();
        newLeaf.taxonId  = x;

        if (start.isLeaf()) {
            // Case 1: stopped at a leaf — graft as sibling.
            // Mirrors ASTRAL-MP's:
            //   newnode = start.getParent().createChild(name);   // append at end
            //   newinternalnode = start.getParent().createChild();  // append at end
            //   newinternalnode.adoptChild(start);   // moves start out of P.children
            //   newinternalnode.adoptChild(newnode); // moves newnode out of P.children
            // Net effect on P.children: start (wherever it was) and newnode are
            // removed, newinternalnode is the only addition at end of list.
            // For our binary representation that means:
            //   P.left  = otherChild (the child that wasn't start)
            //   P.right = newInternal
            // — regardless of which slot start was in.
            TreeNode newInternal = new TreeNode();
            TreeNode p           = start.parent;

            newInternal.left   = start;
            newInternal.right  = newLeaf;
            newInternal.parent = p;
            start.parent       = newInternal;
            newLeaf.parent     = newInternal;

            // p should never be null here: after rerooting the tree has at least
            // anchor on one side and start on the other, so start is never the root.
            if (p != null) {
                if (p.left == start) {
                    // start was at .left → move otherChild (p.right) to .left,
                    // place newInternal at .right (end of list).
                    p.left  = p.right;
                    p.right = newInternal;
                } else {
                    // start was at .right → otherChild already at .left, just replace start.
                    p.right = newInternal;
                }
            }
        } else {
            // Case 2: stopped at an internal node — push children down.
            // c1 and c2 are start.left / start.right as captured by the loop.
            TreeNode newInternal = new TreeNode();

            newInternal.left   = c1;
            newInternal.right  = c2;
            newInternal.parent = start;
            c1.parent          = newInternal;
            c2.parent          = newInternal;

            start.left  = newLeaf;
            start.right = newInternal;
            newLeaf.parent = start;
        }

        return newLeaf;
    }

    /** Insert a taxon while retaining every unresolved child of a polytomy. */
    private static TreeNode insertTaxonPolytomy(int x, TreeNode start,
                                                TreeNode c1, TreeNode c2) {
        TreeNode newLeaf = new TreeNode();
        newLeaf.taxonId = x;

        if (start.isLeaf()) {
            TreeNode parent = start.parent;
            TreeNode newInternal = new TreeNode();
            setBinaryChildren(newInternal, start, newLeaf);
            start.parent = newInternal;
            newLeaf.parent = newInternal;
            removeAndAppend(parent, start, newInternal);
            newInternal.parent = parent;
        } else {
            TreeNode newInternal = new TreeNode();
            setBinaryChildren(newInternal, c1, c2);
            c1.parent = newInternal;
            c2.parent = newInternal;

            if (start.isPolytomous()) {
                TreeNode[] children = start.children;
                int out = 0;
                for (TreeNode child : children) {
                    if (child != c1 && child != c2) children[out++] = child;
                }
                children[out++] = newLeaf;
                children[out++] = newInternal;
                if (out != children.length) {
                    throw new IllegalStateException(
                        "Completion insertion could not locate both navigation children");
                }
                setChildren(start, children);
            } else {
                setBinaryChildren(start, newLeaf, newInternal);
            }
            newLeaf.parent = start;
            newInternal.parent = start;
        }
        return newLeaf;
    }

    // ── Preprocessing reroot ─────────────────────────────────────────────────

    /**
     * Replicate ASTRAL-MP's WQDataCollection.reroot() preprocessing:
     *
     *   1. Find the internal node whose subtree leaf count is closest to
     *      leafCount/2  (the most balanced split).
     *      Mirrors: Math.abs(n - node.getLeafCount()) < dist  with
     *               n = leafCount/2  and the exact signed-dist update so that
     *               tie-breaking matches ASTRAL-MP's post-traversal order.
     *
     *   2. Move every direct leaf-child to the RIGHT slot (= end of list in
     *      ASTRAL-MP's n-ary tree).  Effect: leftmostTaxon() always descends
     *      into an internal child first, matching getLeftmostLeaf() in ASTRAL-MP
     *      after the leaf-to-end reordering.
     *
     *   3. Reroot at the balanced node's edge (if it is not already the root),
     *      using the same path-reversal as rerootAtLeafEdge() — which works for
     *      any node, leaf or internal.
     *
     * Called once per tree before the insertion loop.
     * O(n) time, O(depth) stack space for the recursive leaf-move pass.
     *
     * HINT for n-ary: step 2 should move ALL leaf children to the tail of
     * node.children rather than just swapping left/right.
     */
    private static TreeNode preprocessReroot(TreeNode root, int leafCount) {
        // Step 1: find the most balanced internal node in POST-ORDER.
        // Mirror ASTRAL-MP reroot() exactly:
        //   n    = leafCount / 2
        //   dist = n  (initial, positive)
        //   update only when Math.abs(n - sub) < dist (strict)
        //   dist = n - sub  (signed! — once a node with sub > n is found dist
        //                    goes negative and no further updates happen)
        // The signed update means: the algorithm picks the LAST post-order node
        // that strictly improves |n - sub|, and stops as soon as sub > n is
        // accepted (dist negative ⟹ no |...| can be < negative).
        int half    = leafCount / 2;
        int[] dist  = { half };            // mutable: mirrors ASTRAL-MP's dist variable
        TreeNode[] bestNode = { null };    // mirrors newroot (null = use existing root)

        // Post-order DFS (left, right, node) — mirrors tr.postTraverse()
        findBalancedNodePostOrder(root, half, dist, bestNode);

        // Step 2: move direct leaf-children to the right slot at every internal node.
        // Applied to the ORIGINAL tree (before rerooting), same order as ASTRAL-MP:
        // reroot() collects leaf children and moves them to end in the same traversal
        // that finds the balanced node, then reroots last.
        moveLeafChildToRight(root);

        // Step 3: reroot at best node's edge (if not already the root).
        if (bestNode[0] != null && bestNode[0] != root) {
            root = rerootAtLeafEdge(bestNode[0], root);
            // After: root.left = bestNode, root.right = reversed-path rest-of-tree.
        }

        return root;
    }

    /** ASTRAL-MP-compatible preprocessing over canonical n-ary child lists. */
    private static TreeNode preprocessRerootPolytomy(TreeNode root, int leafCount) {
        int half = leafCount / 2;
        int[] dist = {half};
        TreeNode[] bestNode = {null};
        findBalancedNodePostOrderPolytomy(root, half, dist, bestNode);
        moveLeafChildToEndPolytomy(root);
        if (bestNode[0] != null && bestNode[0] != root) {
            root = rerootAtEdgePolytomy(bestNode[0], root);
        }
        return root;
    }

    private static void findBalancedNodePostOrderPolytomy(TreeNode node, int half,
                                                           int[] dist, TreeNode[] best) {
        if (node.isLeaf()) return;
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) {
                findBalancedNodePostOrderPolytomy(child, half, dist, best);
            }
        } else {
            findBalancedNodePostOrderPolytomy(node.left, half, dist, best);
            findBalancedNodePostOrderPolytomy(node.right, half, dist, best);
        }
        int sub = node.rangeEnd - node.rangeStart;
        if (Math.abs(half - sub) < dist[0]) {
            best[0] = node;
            dist[0] = half - sub;
        }
    }

    private static void moveLeafChildToEndPolytomy(TreeNode node) {
        if (node == null || node.isLeaf()) return;
        if (!node.isPolytomous()) {
            moveLeafChildToEndPolytomy(node.left);
            moveLeafChildToEndPolytomy(node.right);
            if (node.left.isLeaf()) {
                TreeNode tmp = node.left;
                node.left = node.right;
                node.right = tmp;
            }
            return;
        }
        TreeNode[] children = node.children;
        for (TreeNode child : children) moveLeafChildToEndPolytomy(child);
        for (int i = 0; i < children.length; i++) {
            if (children[i].isLeaf()) {
                if (i != children.length - 1) {
                    TreeNode leaf = children[i];
                    System.arraycopy(children, i + 1, children, i, children.length - i - 1);
                    children[children.length - 1] = leaf;
                    setChildren(node, children);
                }
                break;
            }
        }
    }

    /**
     * Post-order DFS to find the most balanced internal node, mirroring
     * ASTRAL-MP's postTraverse() loop in reroot().
     * Updates dist[0] and bestNode[0] in-place.
     */
    private static void findBalancedNodePostOrder(TreeNode node, int half,
                                                   int[] dist, TreeNode[] best) {
        if (node.isLeaf()) return;
        findBalancedNodePostOrder(node.left,  half, dist, best);
        findBalancedNodePostOrder(node.right, half, dist, best);
        // Process this node (post-order = after children)
        int sub = node.rangeEnd - node.rangeStart;
        if (Math.abs(half - sub) < dist[0]) {
            best[0]  = node;
            dist[0]  = half - sub;   // signed — mirrors ASTRAL-MP: dist = n - node.getLeafCount()
        }
    }

    /**
     * Post-order DFS: for each internal node, if left is a leaf and right is
     * not, swap them so the internal child is always on the left.
     *
     * Mirrors ASTRAL-MP's removeChild + createChild(leaf) loop in reroot(),
     * which places every leaf child at the END (= right) of the children list.
     *
     * HINT for n-ary: collect all leaf children, remove and re-append them.
     */
    private static void moveLeafChildToRight(TreeNode node) {
        if (node == null || node.isLeaf()) return;
        moveLeafChildToRight(node.left);
        moveLeafChildToRight(node.right);
        // ASTRAL-MP's reroot() iterates each internal node's children, finds the
        // FIRST leaf child, breaks, and later moves that leaf to the END of the
        // children list (via removeChild + createChild). For our binary tree:
        //   [leaf, leaf]    → [right, left_copy]   (left was first leaf found)
        //   [leaf, internal] → [internal, leaf_copy]
        //   [internal, leaf] → unchanged (leaf already last)
        //   [internal, internal] → unchanged
        // The single rule that reproduces all four cases is:
        //   if the LEFT child is a leaf, swap left and right.
        if (node.left.isLeaf()) {
            TreeNode tmp  = node.left;
            node.left     = node.right;
            node.right    = tmp;
            // parents already point to `node` — no reparent needed since both
            // children were already direct children of `node`.
        }
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    /**
     * Walk left-child pointers until a leaf; return its taxonId.
     *
     * Used both for computing c1rep/c2rep during navigation and as a utility
     * inside deepCopyNodes.
     *
     * Complexity: O(depth).
     *
     * HINT for n-ary: walk children.get(0) instead of node.left.
     */
    private static int leftmostTaxon(TreeNode node) {
        while (!node.isLeaf()) node = node.left;
        return node.taxonId;
    }

    // ── Deep copy ────────────────────────────────────────────────────────────

    /**
     * Recursively deep-copy a TreeNode subtree, also populating taxonNode[].
     *
     * The copies are fresh objects with the same taxonId/rangeStart/rangeEnd
     * but independent parent/left/right pointers — mutations to the copy do
     * not affect the original Tree's nodes (preserving originalTrees' range
     * fields used by PrefixHashArrays and PartitionTable).
     *
     * taxonNode[id] is set to the leaf copy for each leaf encountered.
     *
     * HINT for n-ary: copy node.children list here instead of left/right.
     *
     * @param src       source node
     * @param parent    parent of the copy (null for root)
     * @param taxonNode map from taxon ID → leaf node; populated for leaves
     * @return root of the copied subtree
     */
    private static TreeNode deepCopyNodes(TreeNode src, TreeNode parent,
                                          TreeNode[] taxonNode) {
        if (src == null) return null;
        TreeNode copy    = new TreeNode();
        copy.taxonId     = src.taxonId;
        copy.rangeStart  = src.rangeStart;
        copy.rangeEnd    = src.rangeEnd;
        copy.parent      = parent;
        copy.left        = deepCopyNodes(src.left,  copy, taxonNode);
        copy.right       = deepCopyNodes(src.right, copy, taxonNode);
        if (copy.isLeaf()) {
            taxonNode[copy.taxonId] = copy;
        }
        return copy;
    }

    private static TreeNode deepCopyNodesPolytomy(TreeNode src, TreeNode parent,
                                                   TreeNode[] taxonNode) {
        if (src == null) return null;
        TreeNode copy = new TreeNode();
        copy.taxonId = src.taxonId;
        copy.rangeStart = src.rangeStart;
        copy.rangeEnd = src.rangeEnd;
        copy.parent = parent;
        if (src.isLeaf()) {
            taxonNode[copy.taxonId] = copy;
        } else if (src.isPolytomous()) {
            TreeNode[] children = new TreeNode[src.children.length];
            for (int i = 0; i < children.length; i++) {
                children[i] = deepCopyNodesPolytomy(src.children[i], copy, taxonNode);
            }
            setChildren(copy, children);
        } else {
            copy.left = deepCopyNodesPolytomy(src.left, copy, taxonNode);
            copy.right = deepCopyNodesPolytomy(src.right, copy, taxonNode);
        }
        return copy;
    }

    // ── Rebuild Tree ──────────────────────────────────────────────────────────

    /** Reconstruct a Tree object from the mutated TreeNode structure. */
    private static Tree rebuildTree(int treeIndex, TreeNode root, int n) {
        int[] postorderArray = new int[n];
        int[] counter        = {0};
        assignRangesAndFill(root, postorderArray, counter);
        int leafCount = counter[0];
        postorderArray = Arrays.copyOf(postorderArray, leafCount);

        int[] positionMap = new int[n];
        Arrays.fill(positionMap, -1);
        for (int j = 0; j < leafCount; j++) positionMap[postorderArray[j]] = j;

        return new Tree(treeIndex, root, postorderArray, positionMap, leafCount, n);
    }

    private static Tree rebuildTreePolytomy(int treeIndex, TreeNode root, int n) {
        int[] postorderArray = new int[n];
        int[] counter = {0};
        assignRangesAndFillPolytomy(root, postorderArray, counter);
        int leafCount = counter[0];
        postorderArray = Arrays.copyOf(postorderArray, leafCount);
        int[] positionMap = new int[n];
        Arrays.fill(positionMap, -1);
        for (int j = 0; j < leafCount; j++) positionMap[postorderArray[j]] = j;
        return new Tree(treeIndex, root, postorderArray, positionMap, leafCount, n, true);
    }

    private static void assignRangesAndFill(TreeNode node, int[] arr, int[] counter) {
        if (node.isLeaf()) {
            node.rangeStart  = counter[0];
            node.rangeEnd    = counter[0] + 1;
            arr[counter[0]]  = node.taxonId;
            counter[0]++;
            return;
        }
        assignRangesAndFill(node.left,  arr, counter);
        assignRangesAndFill(node.right, arr, counter);
        node.rangeStart = node.left.rangeStart;
        node.rangeEnd   = node.right.rangeEnd;
    }

    private static void assignRangesAndFillPolytomy(TreeNode node, int[] arr, int[] counter) {
        if (node.isLeaf()) {
            node.rangeStart = counter[0];
            node.rangeEnd = counter[0] + 1;
            arr[counter[0]++] = node.taxonId;
            return;
        }
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) {
                assignRangesAndFillPolytomy(child, arr, counter);
            }
        } else {
            assignRangesAndFillPolytomy(node.left, arr, counter);
            assignRangesAndFillPolytomy(node.right, arr, counter);
        }
        node.rangeStart = node.left.rangeStart;
        node.rangeEnd = node.right.rangeEnd;
    }

    private static TreeNode childAt(TreeNode node, int index) {
        if (node.isPolytomous()) return node.children[index];
        return index == 0 ? node.left : node.right;
    }

    private static void setChildren(TreeNode node, TreeNode[] children) {
        if (children.length < 2) {
            throw new IllegalArgumentException("Internal node must have at least two children");
        }
        node.left = children[0];
        node.right = children[children.length - 1];
        node.children = children.length > 2 ? children : null;
    }

    private static void setBinaryChildren(TreeNode node, TreeNode left, TreeNode right) {
        node.left = left;
        node.right = right;
        node.children = null;
    }

    private static void removeAndAppend(TreeNode node, TreeNode remove, TreeNode append) {
        if (!node.isPolytomous()) {
            TreeNode other;
            if (node.left == remove) other = node.right;
            else if (node.right == remove) other = node.left;
            else throw new IllegalStateException(
                "Completion child-list mutation lost its path child");
            node.left = other;
            node.right = append;
            return;
        }
        TreeNode[] children = node.children;
        int out = 0;
        boolean found = false;
        for (TreeNode child : children) {
            if (child == remove) {
                found = true;
            } else {
                children[out++] = child;
            }
        }
        if (!found || out != children.length - 1) {
            throw new IllegalStateException("Completion child-list mutation lost its path child");
        }
        children[out] = append;
        setChildren(node, children);
    }

    private static void replaceChild(TreeNode node, TreeNode remove, TreeNode replacement) {
        if (!node.isPolytomous()) {
            if (node.left == remove) node.left = replacement;
            else if (node.right == remove) node.right = replacement;
            else throw new IllegalStateException(
                "Completion child-list mutation lost the old root");
            return;
        }
        TreeNode[] children = node.children;
        for (int i = 0; i < children.length; i++) {
            if (children[i] == remove) {
                children[i] = replacement;
                setChildren(node, children);
                return;
            }
        }
        throw new IllegalStateException("Completion child-list mutation lost the old root");
    }

    private static void requireBinaryRoot(TreeNode root) {
        if (root.isPolytomous()) {
            throw new IllegalStateException("Completion reroot expects a binary root, found degree "
                + root.children.length);
        }
    }
}
