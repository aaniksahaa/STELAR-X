package stelarx.greedy;

import java.util.ArrayList;
import java.util.List;

/**
 * Incremental laminar forest used by greedy consensus construction.
 *
 * State at any moment is a rooted tree (a "partial consensus tree"):
 *   - virtualRoot is a single node containing all n taxa as descendants
 *   - initially virtualRoot has the n leaves as direct children (star tree)
 *   - each accepted INSERT(C) creates a new internal node that adopts a
 *     sub-collection of an existing node's children
 *
 * The forest exposes:
 *   - the {@link LaminarNode} graph (parent, children, size)
 *   - {@link #leafNode} for taxon-id → leaf node lookup
 *   - {@link #lcaOfTaxa}: deepest common ancestor of a set of taxa
 *     (matches what ASTRAL-MP's SchieberVishkinLCA.getLCA computes)
 *   - {@link #childOnPathTo}: the child of v that is an ancestor of taxon t
 *
 * LCA-of-set is computed with epoch-based path marking — no node depths
 * are maintained (which would otherwise need O(|subtree|) updates per
 * INSERT, since adopting a child shifts that subtree one level deeper).
 *
 * Capacity grows on demand: the max possible node count for a laminar
 * forest over n taxa is n leaves + 1 virtual root + (n − 1) internal
 * absorbers = 2n.  We pre-size the parallel mark arrays to 2n+1.
 */
final class LaminarForest {

    final int numTaxa;
    final LaminarNode virtualRoot;
    /** taxonId → leaf node. */
    private final LaminarNode[] leafNodes;
    /** All allocated nodes in id order. */
    private final List<LaminarNode> nodes;

    // ── Epoch-based LCA scratch storage (sized to 2n+1, grown if exceeded) ──
    private int[] markEpoch;      // markEpoch[id] == currentEpoch ⇒ node is on chain1
    private int[] markPos;        // when marked, position in chain1 (0 = leaf-side)
    private int   epoch = 0;

    LaminarForest(int numTaxa) {
        if (numTaxa < 1) throw new IllegalArgumentException("numTaxa < 1");
        this.numTaxa = numTaxa;
        this.nodes = new ArrayList<>(2 * numTaxa + 1);

        // Virtual root: parent = -1, size = numTaxa, internal
        this.virtualRoot = new LaminarNode(0, -1, numTaxa, -1);
        this.nodes.add(virtualRoot);

        // One leaf node per taxon, direct child of virtualRoot
        this.leafNodes = new LaminarNode[numTaxa];
        for (int t = 0; t < numTaxa; t++) {
            LaminarNode leaf = new LaminarNode(/*id*/ t + 1, /*parent*/ 0,
                                               /*size*/ 1, /*taxon*/ t);
            this.nodes.add(leaf);
            this.leafNodes[t] = leaf;
            this.virtualRoot.children.add(leaf);
        }

        int cap = 2 * numTaxa + 4;
        this.markEpoch = new int[cap];
        this.markPos   = new int[cap];
    }

    int nodeCount() { return nodes.size(); }
    LaminarNode leafNode(int taxonId) { return leafNodes[taxonId]; }
    LaminarNode getNode(int id)       { return nodes.get(id); }

    /**
     * Allocate a new internal node, attach as a child of {@code parent}.
     * Size and child wiring are the caller's responsibility (the caller is
     * INSERT, which moves a subset of parent's children to the new node).
     */
    LaminarNode createInternal(LaminarNode parent) {
        int id = nodes.size();
        LaminarNode n = new LaminarNode(id, parent.id, 0, -1);
        nodes.add(n);
        ensureMarkCapacity(id + 1);
        return n;
    }

    private void ensureMarkCapacity(int needed) {
        if (needed <= markEpoch.length) return;
        int cap = Math.max(needed, markEpoch.length * 2);
        int[] e2 = new int[cap];
        int[] p2 = new int[cap];
        System.arraycopy(markEpoch, 0, e2, 0, markEpoch.length);
        System.arraycopy(markPos,   0, p2, 0, markPos.length);
        this.markEpoch = e2;
        this.markPos   = p2;
    }

    /**
     * Deepest common ancestor of the given taxa in the laminar forest.
     * Equivalent to ASTRAL-MP's {@code SchieberVishkinLCA.getLCA(leafNodes)}.
     * Time: O(|taxa| × tree depth) amortized.
     *
     * Pre: taxa.length >= 1; all entries are valid taxon ids.
     */
    LaminarNode lcaOfTaxa(int[] taxa) {
        if (taxa.length == 0) return virtualRoot;
        epoch++;

        // ── Pass 1: mark the leaf→root chain for taxa[0] with positions ──
        // chain1[k] is the k-th node from leaf upward; pos 0 = the leaf itself
        List<LaminarNode> chain1 = new ArrayList<>(8);
        LaminarNode cur = leafNodes[taxa[0]];
        int pos = 0;
        while (cur != null) {
            int id = cur.id;
            markEpoch[id] = epoch;
            markPos[id]   = pos++;
            chain1.add(cur);
            cur = (cur.parent < 0) ? null : nodes.get(cur.parent);
        }

        if (taxa.length == 1) {
            // LCA of a single leaf = that leaf
            return leafNodes[taxa[0]];
        }

        // ── Pass 2: for each subsequent leaf, walk up to first marked ancestor
        //            (== pairwise LCA with taxa[0]); track the SHALLOWEST one
        //            (largest pos in chain1). ──
        int shallowestPos = -1;
        for (int i = 1; i < taxa.length; i++) {
            LaminarNode c = leafNodes[taxa[i]];
            while (markEpoch[c.id] != epoch) {
                int p = c.parent;
                if (p < 0) {
                    // Should be unreachable: virtual root was marked in pass 1
                    return virtualRoot;
                }
                c = nodes.get(p);
            }
            int p = markPos[c.id];
            if (p > shallowestPos) shallowestPos = p;
        }

        return chain1.get(shallowestPos);
    }

    /**
     * Return the child of {@code v} that is an ancestor (or self) of taxon
     * {@code t}.  Pre: t is in v's subtree (v is an ancestor of leaf t).
     *
     * Walks up from leaf t along parent pointers and returns the first node
     * whose parent is v.  Time: O(distance from leaf to v's child) per call.
     */
    LaminarNode childOnPathTo(LaminarNode v, int taxonId) {
        LaminarNode cur = leafNodes[taxonId];
        while (cur.parent != v.id) {
            int p = cur.parent;
            if (p < 0) {
                // taxon not under v — caller violated precondition
                throw new IllegalStateException(
                    "childOnPathTo: taxon " + taxonId + " not under node " + v.id);
            }
            cur = nodes.get(p);
        }
        return cur;
    }

    /**
     * Adopt the listed children into a new internal node {@code newParent}
     * (which must already be attached as a child of {@code oldParent}, with
     * an empty children list).  Updates parent pointers, the children list
     * on both sides, and the size field of newParent.
     */
    void moveChildren(LaminarNode oldParent, LaminarNode newParent,
                      List<LaminarNode> toMove) {
        // Remove from oldParent in a single linear pass
        // (toMove may be much smaller than oldParent.children, so O(|oldP.children|)
        //  per insert is fine — net O(B̄) across the build because we only ever
        //  re-touch each child at most once per insert where it is moved).
        if (toMove.size() > 0) {
            // Build a fast membership predicate
            java.util.HashSet<LaminarNode> set =
                new java.util.HashSet<>(toMove.size() * 2);
            for (LaminarNode c : toMove) set.add(c);
            List<LaminarNode> kept = new ArrayList<>(oldParent.children.size());
            for (LaminarNode c : oldParent.children) {
                if (!set.contains(c)) kept.add(c);
            }
            oldParent.children.clear();
            oldParent.children.addAll(kept);
        }

        int totalSize = 0;
        for (LaminarNode c : toMove) {
            c.parent = newParent.id;
            newParent.children.add(c);
            totalSize += c.size;
        }
        newParent.size = totalSize;
    }
}
