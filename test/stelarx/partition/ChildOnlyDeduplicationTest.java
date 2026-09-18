package stelarx.partition;

import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** Structural regression for child-only deduplication on incomplete rooted trees. */
public final class ChildOnlyDeduplicationTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass the dedup fixture");

        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        try {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> trees = TreeParser.parseGeneTrees(args[0], registry, true);
            TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 0x5E1A7L);
            PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);
            PartitionTable table = new PartitionTable(trees, pref);

            check(table.size() == 6,
                "expected 6 unique contributing child partitions, observed " + table.size());

            int totalFrequency = 0;
            int polytomyEntries = 0;
            List<Integer> frequencies = new ArrayList<>();
            for (PartitionTable.Entry entry : table.entries()) {
                totalFrequency += entry.frequency;
                frequencies.add(entry.frequency);
                if (entry.exemplar.isPolytomous()) polytomyEntries++;
                check(!(entry.exemplar.d == 3
                        && entry.exemplar.size1 == 1 && entry.exemplar.size2 == 1),
                    "zero-contribution binary cherry reached the table");
            }
            Collections.sort(frequencies);
            check(totalFrequency == 8,
                "expected 8 contributing node occurrences, observed " + totalFrequency);
            check(frequencies.equals(List.of(1, 1, 1, 1, 2, 2)),
                "child partitions under different complements did not merge: " + frequencies);
            check(polytomyEntries == 1,
                "same polytomy children under different complements did not merge");

            System.out.println("Child-only incomplete-tree deduplication: PASS "
                + "(14 internal nodes -> 8 contributing occurrences -> 6 unique)");
        } finally {
            Threading.shutdown();
        }
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
