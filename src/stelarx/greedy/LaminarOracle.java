package stelarx.greedy;

import stelarx.cluster.Cluster;
import stelarx.tree.Tree;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * Reference implementation of laminar-forest INSERT — the "obviously correct"
 * brute-force algorithm used to validate {@link LaminarBuilder}.
 *
 * Maintains every node's leaf-id set explicitly (HashSet&lt;Integer&gt;).  LCA is
 * found by walking the tree top-down, descending into whichever child wholly
 * contains C, repeating until no further descent is possible.  Whole-node test
 * is then standard set inclusion.
 *
 * Time per INSERT: O(n × tree_size) — far slower than {@link LaminarBuilder},
 * but transparently correct.  Used only by the verifier.
 */
final class LaminarOracle {

    /** Tree node with explicit leaf set. */
    static final class Node {
        final int id;
        Node parent;
        final Set<Integer> leaves;
        final List<Node> children = new ArrayList<>();
        final int taxonId;  // -1 if internal

        Node(int id, Node parent, Set<Integer> leaves, int taxonId) {
            this.id = id; this.parent = parent;
            this.leaves = leaves; this.taxonId = taxonId;
        }

        boolean isLeaf() { return taxonId >= 0; }
    }

    private final int numTaxa;
    private final List<Tree> trees;
    private final Node virtualRoot;
    private int nextId = 0;

    LaminarOracle(int numTaxa, List<Tree> trees) {
        this.numTaxa = numTaxa;
        this.trees = trees;
        Set<Integer> all = new HashSet<>();
        for (int t = 0; t < numTaxa; t++) all.add(t);
        this.virtualRoot = new Node(nextId++, null, all, -1);
        for (int t = 0; t < numTaxa; t++) {
            Set<Integer> single = new HashSet<>();
            single.add(t);
            Node leaf = new Node(nextId++, virtualRoot, single, t);
            virtualRoot.children.add(leaf);
        }
    }

    Node virtualRoot() { return virtualRoot; }

    /** Brute-force insert: returns the same Outcome enum as LaminarBuilder. */
    LaminarBuilder.Outcome insert(Bipartition b) {
        if (b.size <= 1 || b.size == numTaxa) return LaminarBuilder.Outcome.SKIP_TRIVIAL;
        Set<Integer> C = enumerateTaxaSet(b.canonicalExemplar);

        // Find LCA: deepest node whose leaves ⊇ C.
        Node lca = findLCA(C);
        if (lca == null) {
            // C is not a subset of any node — impossible in a tree with virtualRoot
            // covering all taxa.  Defensive.
            return LaminarBuilder.Outcome.REJECT_CROSS_CUT;
        }

        // Test each child of LCA: child wholly inside C, or disjoint, or split?
        List<Node> moved = new ArrayList<>();
        int accountedFor = 0;
        for (Node child : lca.children) {
            int overlap = 0;
            for (int t : C) {
                if (child.leaves.contains(t)) overlap++;
            }
            if (overlap == 0) continue;                         // disjoint
            if (overlap != child.leaves.size()) {               // split
                return LaminarBuilder.Outcome.REJECT_CROSS_CUT;
            }
            moved.add(child);
            accountedFor += overlap;
        }
        if (accountedFor != C.size()) {
            return LaminarBuilder.Outcome.REJECT_CROSS_CUT;
        }
        if (moved.size() < 2 || moved.size() == lca.children.size()) {
            return LaminarBuilder.Outcome.SKIP_REDUNDANT;
        }

        // Create new node, adopt moved children.
        Set<Integer> wLeaves = new HashSet<>();
        for (Node c : moved) wLeaves.addAll(c.leaves);
        Node w = new Node(nextId++, lca, wLeaves, -1);
        lca.children.removeAll(moved);
        lca.children.add(w);
        for (Node c : moved) {
            c.parent = w;
            w.children.add(c);
        }
        return LaminarBuilder.Outcome.ACCEPT;
    }

    /** Deepest node whose leaf set ⊇ C.  Top-down descent. */
    private Node findLCA(Set<Integer> C) {
        Node cur = virtualRoot;
        outer:
        while (true) {
            for (Node child : cur.children) {
                if (child.leaves.containsAll(C)) {
                    cur = child;
                    continue outer;
                }
            }
            return cur;
        }
    }

    private Set<Integer> enumerateTaxaSet(Cluster ex) {
        Tree tree = trees.get(ex.treeIndex);
        int[] arr = tree.postorderArray;
        Set<Integer> out = new HashSet<>(ex.size * 2);
        if (!ex.complement) {
            for (int i = ex.left; i < ex.right; i++) out.add(arr[i]);
        } else {
            // Super-complement w.r.t. ALL n taxa — must include taxa missing
            // from the source tree (positionMap[t] < 0).
            int[] pos = tree.positionMap;
            for (int t = 0; t < numTaxa; t++) {
                int p = pos[t];
                if (p < 0 || p < ex.left || p >= ex.right) out.add(t);
            }
        }
        return out;
    }

    // ── Snapshot for cross-checking ─────────────────────────────────────────

    /**
     * Emit a Newick-style canonical string of the current tree.
     * Children are sorted by their sorted-leaf-id-tuple so two structurally
     * identical trees produce identical strings regardless of insertion order.
     */
    String canonicalNewick() {
        return buildCanonicalNewick(virtualRoot) + ";";
    }

    private String buildCanonicalNewick(Node n) {
        if (n.isLeaf()) return Integer.toString(n.taxonId);
        // Sort children by their min-leaf-id, then by leaf-count, then by sorted-leaf-list
        List<Node> sortedChildren = new ArrayList<>(n.children);
        sortedChildren.sort((a, bb) -> {
            int la = minLeaf(a), lb = minLeaf(bb);
            if (la != lb) return Integer.compare(la, lb);
            int sa = a.leaves.size(), sb = bb.leaves.size();
            return Integer.compare(sa, sb);
        });
        StringBuilder sb = new StringBuilder("(");
        for (int i = 0; i < sortedChildren.size(); i++) {
            if (i > 0) sb.append(',');
            sb.append(buildCanonicalNewick(sortedChildren.get(i)));
        }
        sb.append(')');
        return sb.toString();
    }

    private int minLeaf(Node n) {
        int m = Integer.MAX_VALUE;
        for (int t : n.leaves) if (t < m) m = t;
        return m;
    }

    /** Sorted leaf-set string representation per node (for set-equality checks). */
    String canonicalLeafSets() {
        List<String> lines = new ArrayList<>();
        collect(virtualRoot, lines);
        java.util.Collections.sort(lines);
        return String.join("\n", lines);
    }

    private void collect(Node n, List<String> out) {
        if (!n.isLeaf() && n != virtualRoot) {
            TreeSet<Integer> sorted = new TreeSet<>(n.leaves);
            StringBuilder sb = new StringBuilder("{");
            boolean first = true;
            for (int t : sorted) {
                if (!first) sb.append(',');
                sb.append(t);
                first = false;
            }
            sb.append('}');
            out.add(sb.toString());
        }
        for (Node c : n.children) collect(c, out);
    }
}
