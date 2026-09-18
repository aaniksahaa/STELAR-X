package stelarx.cluster;

import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.ArrayList;
import java.util.List;

/** Differential check for the allocation-free cross-tree residual probe. */
public final class ResidualLookupTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass one gene-tree fixture");

        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        try {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> trees = TreeParser.parseGeneTrees(args[0], registry, false);
            TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 1L);
            PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);
            ClusterTable table = new ClusterTable(trees, pref, registry.size());

            List<ClusterHash> parents = new ArrayList<>();
            parents.add(table.getAllTaxaHash());
            for (ClusterTable.Entry entry : table.entries()) parents.add(entry.hash);

            long checks = 0;
            for (ClusterHash parent : parents) {
                ClusterTable.ResidualLookup lookup = table.newResidualLookup();
                for (ClusterTable.Entry childEntry : table.entries()) {
                    ClusterHash child = childEntry.hash;
                    ClusterTable.Entry expected = table.get(ClusterHash.residual(parent, child));
                    ClusterHash actual = lookup.find(parent, child);
                    if ((expected == null) != (actual == null)
                            || (expected != null && !expected.hash.equals(actual))) {
                        throw new AssertionError("residual lookup mismatch for "
                            + parent + " minus " + child);
                    }
                    checks++;
                }
            }
            System.out.println("Allocation-free residual lookup: PASS (" + checks
                + " exact comparisons)");
        } finally {
            Threading.shutdown();
        }
    }
}
