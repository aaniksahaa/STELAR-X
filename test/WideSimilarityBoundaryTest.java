import stelarx.completion.EulerTourBuilder;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.Random;

/** Boundary and widened-child-size checks without allocating an n-by-n matrix. */
public final class WideSimilarityBoundaryTest {
    private static TreeNode buildBalanced(int lo, int hi, int[] postorder) {
        if (hi - lo == 1) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId = lo;
            leaf.rangeStart = lo;
            leaf.rangeEnd = hi;
            postorder[lo] = lo;
            return leaf;
        }
        int mid = (lo + hi) >>> 1;
        TreeNode node = new TreeNode();
        node.left = buildBalanced(lo, mid, postorder);
        node.right = buildBalanced(mid, hi, postorder);
        node.left.parent = node;
        node.right.parent = node;
        node.rangeStart = lo;
        node.rangeEnd = hi;
        return node;
    }

    private static Tree buildTree(int n, int leftRootSize) {
        int[] postorder = new int[n];
        int[] position = new int[n];
        for (int i = 0; i < n; i++) position[i] = i;
        TreeNode root;
        if (leftRootSize > 0) {
            root = new TreeNode();
            root.left = buildBalanced(0, leftRootSize, postorder);
            root.right = buildBalanced(leftRootSize, n, postorder);
            root.left.parent = root;
            root.right.parent = root;
            root.rangeStart = 0;
            root.rangeEnd = n;
        } else {
            root = buildBalanced(0, n, postorder);
        }
        return new Tree(0, root, postorder, position, n, n);
    }

    private static int directArgmin(int[] depths, int lo, int hi) {
        int best = lo;
        for (int p = lo + 1; p <= hi; p++) {
            if (depths[p] < depths[best]) best = p;
        }
        return best;
    }

    private static void checkRange(EulerTourBuilder.WideTourData wide, int lo, int hi) {
        int expected = directArgmin(wide.depths, lo, hi);
        int actual = EulerTourBuilder.queryWideArgmin(wide, lo, hi);
        if (actual != expected) {
            throw new AssertionError("wide boundary RMQ mismatch for [" + lo + "," + hi
                + "]: expected=" + expected + " actual=" + actual);
        }
    }

    public static void main(String[] args) {
        Tree atLimit = buildTree(21_846, 0);
        if (EulerTourBuilder.tourLength(atLimit) != 65_536) {
            throw new AssertionError("compact boundary length is not 65,536");
        }
        EulerTourBuilder.buildFull(atLimit, atLimit.leafCount);

        Tree overLimit = buildTree(21_847, 0);
        if (EulerTourBuilder.tourLength(overLimit) != 65_539) {
            throw new AssertionError("wide boundary length is not 65,539");
        }
        try {
            EulerTourBuilder.buildFull(overLimit, overLimit.leafCount);
            throw new AssertionError("compact builder accepted a tour beyond 65,536");
        } catch (IllegalArgumentException expected) {
            if (!expected.getMessage().contains("compact 16-bit RMQ")) throw expected;
        }

        EulerTourBuilder.WideTourData wide =
            EulerTourBuilder.buildWide(overLimit, overLimit.leafCount);
        Random random = new Random(0xA57A1L);
        for (int check = 0; check < 200; check++) {
            int a = random.nextInt(wide.tourLen);
            int b = random.nextInt(wide.tourLen);
            int lo = Math.min(a, b), hi = Math.max(a, b);
            checkRange(wide, lo, hi);
        }
        for (int boundary = EulerTourBuilder.WIDE_BLOCK_SIZE;
                boundary < wide.tourLen;
                boundary += EulerTourBuilder.WIDE_BLOCK_SIZE) {
            checkRange(wide, boundary - 1, boundary);
            checkRange(wide, Math.max(0, boundary - 257),
                Math.min(wide.tourLen - 1, boundary + 257));
        }
        checkRange(wide, 0, wide.tourLen - 1);
        checkRange(wide, wide.tourLen - 173, wide.tourLen - 1);

        Tree largeChild = buildTree(40_000, 33_000);
        EulerTourBuilder.WideTourData childWide =
            EulerTourBuilder.buildWide(largeChild, largeChild.leafCount);
        int rootPos = EulerTourBuilder.queryWideArgmin(childWide,
            childWide.firstOcc[0], childWide.firstOcc[39_999]);
        if (childWide.eulerLeftChildS[rootPos] != 33_000
                || childWide.eulerRightChildS[rootPos] != 7_000) {
            throw new AssertionError("wide child sizes were truncated: left="
                + childWide.eulerLeftChildS[rootPos] + " right="
                + childWide.eulerRightChildS[rootPos]);
        }

        System.out.println("Wide similarity boundaries: PASS");
    }
}
