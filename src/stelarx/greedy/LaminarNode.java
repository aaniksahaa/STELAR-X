package stelarx.greedy;

import java.util.ArrayList;
import java.util.List;

/**
 * One node in a {@link LaminarForest}.
 *
 *   id     — dense integer in [0, forest.nodeCount())
 *   parent — id of laminar-forest parent (-1 for the virtual root)
 *   size   — number of leaves in this node's subtree
 *   taxonId — taxon id for leaves; -1 for internal nodes
 *
 * Children are stored as a mutable list so INSERT can move them between
 * sibling sets in O(degree) per insert.  Total node count is O(n) (≤ 2n)
 * across the whole build, so per-node ArrayList overhead is acceptable.
 */
final class LaminarNode {
    final int id;
    int parent;
    int size;
    final int taxonId;            // -1 if internal
    final List<LaminarNode> children = new ArrayList<>(2);

    LaminarNode(int id, int parent, int size, int taxonId) {
        this.id      = id;
        this.parent  = parent;
        this.size    = size;
        this.taxonId = taxonId;
    }

    boolean isLeaf()        { return taxonId >= 0; }
}
