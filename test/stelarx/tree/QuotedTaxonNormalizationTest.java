package stelarx.tree;

import stelarx.taxon.TaxonRegistry;

import java.util.Arrays;
import java.util.List;

/** Regression: valid Newick quoting must not change a taxon's identity. */
public final class QuotedTaxonNormalizationTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 2) throw new IllegalArgumentException("pass both quote fixtures");

        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(args[0], registry, true);

        check(trees.size() == 2, "expected two trees");
        check(registry.size() == 4,
            "quoted and unquoted spellings registered different taxa: " + registry.size());
        int bird = registry.findId("Acanthisitta_chloris");
        check(bird >= 0, "normalized underscore taxon was not registered");
        check(registry.findId("'Acanthisitta_chloris'") < 0,
            "literal quote characters leaked into the taxon registry");
        check(Arrays.equals(trees.get(0).postorderArray, trees.get(1).postorderArray),
            "equivalent quoted and unquoted trees parsed differently");

        for (Tree tree : trees) {
            String newick = tree.toNewick(registry);
            check(newick.contains("Acanthisitta_chloris"), "taxon missing from output");
            check(!newick.contains("'Acanthisitta_chloris'"),
                "ordinary underscore taxon was unnecessarily quoted: " + newick);
        }

        TaxonRegistry escapedRegistry = new TaxonRegistry();
        List<Tree> escapedTrees = TreeParser.parseGeneTrees(args[1], escapedRegistry, true);
        check(escapedRegistry.findId("O'Brien_species") >= 0,
            "doubled quote in a quoted Newick label was not decoded");
        check(escapedTrees.get(0).toNewick(escapedRegistry).contains("'O''Brien_species'"),
            "special taxon name was not safely re-encoded");

        System.out.println("Quoted/unquoted Java taxon normalization: PASS");
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
