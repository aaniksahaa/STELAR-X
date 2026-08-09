package core;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.Arrays;
import java.util.Random;
import org.junit.jupiter.api.Test;

class InverseIndexManagerTest {

    @Test
    void fusedIntersectionsMatchFourIndependentIntersections() {
        int[][] orderings = {
            {0, 1, 2, 3, 4, 5, 6, 7},
            {7, 2, 5, 0, 6, 3, 1, 4},
            {6, 0, 4, 2, 7} // deliberately incomplete tree
        };
        InverseIndexManager manager = new InverseIndexManager(orderings, 8);

        assertFusedMatchesLegacy(manager, 0, 0, 3, 3, 8, 1, 0, 2, 2, 6);
        assertFusedMatchesLegacy(manager, 0, 1, 2, 2, 4, 2, 0, 2, 2, 5);
        assertFusedMatchesLegacy(manager, 2, 0, 1, 1, 3, 1, 0, 4, 4, 8);
    }

    @Test
    void randomizedFusedIntersectionsMatchNaiveReference() {
        Random random = new Random(0x5e1aL);
        int taxa = 100;
        int[][] orderings = new int[12][];
        for (int tree = 0; tree < orderings.length; tree++) {
            int present = 25 + random.nextInt(76);
            int[] permutation = shuffledTaxa(taxa, random);
            orderings[tree] = Arrays.copyOf(permutation, present);
        }
        InverseIndexManager manager = new InverseIndexManager(orderings, taxa);

        for (int iteration = 0; iteration < 2_000; iteration++) {
            int tree1 = random.nextInt(orderings.length);
            int tree2 = random.nextInt(orderings.length);
            int[] cuts1 = cuts(orderings[tree1].length, random);
            int[] cuts2 = cuts(orderings[tree2].length, random);
            assertFusedMatchesLegacy(manager,
                tree1, cuts1[0], cuts1[1], cuts1[1], cuts1[2],
                tree2, cuts2[0], cuts2[1], cuts2[1], cuts2[2]);
        }
    }

    @Test
    void rejectsDuplicateTaxaInAnOrdering() {
        assertThrows(IllegalArgumentException.class,
            () -> new InverseIndexManager(new int[][] {{0, 1, 1}}, 3));
    }

    private static void assertFusedMatchesLegacy(InverseIndexManager manager,
            int tree1, int ls1, int le1, int rs1, int re1,
            int tree2, int ls2, int le2, int rs2, int re2) {
        int[] expected = {
            manager.getRangeIntersectionSize(tree1, ls1, le1, tree2, ls2, le2),
            manager.getRangeIntersectionSize(tree1, ls1, le1, tree2, rs2, re2),
            manager.getRangeIntersectionSize(tree1, rs1, re1, tree2, ls2, le2),
            manager.getRangeIntersectionSize(tree1, rs1, re1, tree2, rs2, re2)
        };
        int[] actual = new int[4];
        manager.getBipartitionIntersections(
            tree1, ls1, le1, rs1, re1, tree2, ls2, le2, rs2, re2, actual);
        assertArrayEquals(expected, actual);
    }

    private static int[] shuffledTaxa(int taxa, Random random) {
        int[] values = new int[taxa];
        for (int i = 0; i < taxa; i++) values[i] = i;
        for (int i = taxa - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            int tmp = values[i]; values[i] = values[j]; values[j] = tmp;
        }
        return values;
    }

    private static int[] cuts(int length, Random random) {
        int start = random.nextInt(length);
        int middle = start + 1 + random.nextInt(length - start);
        int end = middle + random.nextInt(length - middle + 1);
        return new int[] {start, middle, end};
    }
}
