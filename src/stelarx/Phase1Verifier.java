package stelarx;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.io.*;
import java.util.Arrays;
import java.util.List;

/**
 * Dumps Phase-1 parse results in a structured, human-readable form for
 * manual verification against test/expected/ files.
 */
public class Phase1Verifier {

    public static void dump(List<Tree> trees, TaxonRegistry registry, String outFile)
            throws IOException {
        PrintStream out = (outFile != null)
            ? new PrintStream(new FileOutputStream(outFile))
            : System.out;

        out.printf("=== Phase 1 Parse Verification ===%n");
        out.printf("Total taxa: %d%n", registry.size());
        out.printf("Total trees: %d%n%n", trees.size());

        // Taxon ID mapping
        out.println("--- Taxon ID assignments ---");
        for (int id = 0; id < registry.size(); id++) {
            out.printf("  %s = %d%n", registry.getName(id), id);
        }
        out.println();

        // Per-tree details
        for (Tree t : trees) {
            out.printf("--- Tree %d ----%n", t.treeIndex);
            out.printf("  leafCount  : %d%n", t.leafCount);
            out.printf("  isComplete : %b%n", t.isComplete);
            out.printf("  postorder  : %s%n", taxonNames(t.postorderArray, registry));
            out.printf("  postorder  : %s  (IDs)%n", Arrays.toString(t.postorderArray));
            out.println("  positionMap:");
            for (int id = 0; id < registry.size(); id++) {
                int pos = t.positionMap[id];
                if (pos != -1)
                    out.printf("    %s(id=%d) -> pos %d%n", registry.getName(id), id, pos);
                else
                    out.printf("    %s(id=%d) -> ABSENT%n", registry.getName(id), id);
            }
            out.println("  node ranges (DFS pre-order):");
            printNodeRanges(t.root, "    ", out, registry);
            // Sanity checks
            out.println("  --- assertions ---");
            boolean ok = true;
            // 1. positionMap is inverse of postorderArray
            for (int j = 0; j < t.leafCount; j++) {
                int tid = t.postorderArray[j];
                if (t.positionMap[tid] != j) {
                    out.printf("  FAIL: positionMap[%s] = %d, expected %d%n",
                        registry.getName(tid), t.positionMap[tid], j);
                    ok = false;
                }
            }
            // 2. Root range covers all leaves
            if (t.root.rangeStart != 0 || t.root.rangeEnd != t.leafCount) {
                out.printf("  FAIL: root range [%d,%d) != [0,%d)%n",
                    t.root.rangeStart, t.root.rangeEnd, t.leafCount);
                ok = false;
            }
            // 3. Every internal node: rangeEnd = left.rangeEnd, rangeStart = left.rangeStart
            //    and right.rangeStart = left.rangeEnd (contiguous)
            boolean rangeOk = checkRanges(t.root, out);
            if (!rangeOk) ok = false;

            if (ok) out.println("  ALL ASSERTIONS PASSED");
            out.println();
        }

        if (outFile != null) out.close();
    }

    // -------------------------------------------------------------------------

    private static void printNodeRanges(TreeNode n, String indent,
                                         PrintStream out, TaxonRegistry reg) {
        if (n.isLeaf()) {
            out.printf("%sleaf %s : [%d,%d)%n",
                indent, reg.getName(n.taxonId), n.rangeStart, n.rangeEnd);
        } else {
            out.printf("%sinternal%s [%d,%d)%n", indent,
                n.isPolytomous() ? " degree=" + n.children.length : "",
                n.rangeStart, n.rangeEnd);
            if (n.isPolytomous()) {
                for (TreeNode child : n.children) {
                    printNodeRanges(child, indent + "  ", out, reg);
                }
            } else {
                printNodeRanges(n.left,  indent + "  ", out, reg);
                printNodeRanges(n.right, indent + "  ", out, reg);
            }
        }
    }

    private static boolean checkRanges(TreeNode n, PrintStream out) {
        if (n.isLeaf()) return true;
        boolean ok = true;
        TreeNode[] children = n.isPolytomous()
            ? n.children : new TreeNode[]{n.left, n.right};
        // Children must be contiguous and together cover parent.
        for (int i = 1; i < children.length; i++) {
            if (children[i - 1].rangeEnd != children[i].rangeStart) {
                out.printf("  FAIL: children %d/%d not contiguous: %d != %d%n",
                    i - 1, i, children[i - 1].rangeEnd, children[i].rangeStart);
                ok = false;
            }
        }
        if (n.rangeStart != n.left.rangeStart) {
            out.printf("  FAIL: node.rangeStart=%d != left.rangeStart=%d%n",
                n.rangeStart, n.left.rangeStart);
            ok = false;
        }
        if (n.rangeEnd != n.right.rangeEnd) {
            out.printf("  FAIL: node.rangeEnd=%d != right.rangeEnd=%d%n",
                n.rangeEnd, n.right.rangeEnd);
            ok = false;
        }
        for (TreeNode child : children) ok = checkRanges(child, out) && ok;
        return ok;
    }

    private static String taxonNames(int[] ids, TaxonRegistry reg) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < ids.length; i++) {
            if (i > 0) sb.append(", ");
            sb.append(reg.getName(ids[i]));
        }
        sb.append("]");
        return sb.toString();
    }
}
