import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.List;

/** Default-refinement parity with clean.py --deterministic and keep-path coverage. */
public final class PolytomyPreprocessingTest {
    private static final String[] CLEAN_PY_DETERMINISTIC = {
        "(F,((C,(A,B)),(D,E)));",
        "((F,G),(A,(D,(B,C))));",
        "((D,G),((E,(A,C)),B));",
        "(F,((A,D),(G,(B,E))));"
    };

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass the polytomy fixture");

        Parsed serial = parse(args[0], false);
        checkExpected(serial, "serial");

        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        Parsed parallel;
        try {
            parallel = parse(args[0], false);
        } finally {
            Threading.shutdown();
        }
        checkExpected(parallel, "parallel");
        for (int i = 0; i < serial.trees.size(); i++) {
            check(serial.newick(i).equals(parallel.newick(i)),
                "parallel preprocessing changed tree " + i);
        }

        Parsed kept = parse(args[0], true);
        for (Tree tree : kept.trees) {
            check(tree.hasPolytomy,
                "polytomy-preserving parse lost native node in tree " + tree.treeIndex);
            check(!tree.isComplete, "keep-path fixture tree must remain incomplete");
        }
        System.out.println("Polytomy preprocessing/clean.py parity: PASS");
    }

    private static Parsed parse(String path, boolean keep) throws Exception {
        TaxonRegistry registry = new TaxonRegistry();
        return new Parsed(TreeParser.parseGeneTrees(path, registry, keep), registry);
    }

    private static void checkExpected(Parsed parsed, String mode) {
        check(parsed.trees.size() == CLEAN_PY_DETERMINISTIC.length, mode + " tree count");
        for (int i = 0; i < parsed.trees.size(); i++) {
            Tree tree = parsed.trees.get(i);
            check(!tree.hasPolytomy, mode + " tree " + i + " was not fully refined");
            check(CLEAN_PY_DETERMINISTIC[i].equals(parsed.newick(i)),
                mode + " clean.py mismatch for tree " + i + ": " + parsed.newick(i));
        }
    }

    private record Parsed(List<Tree> trees, TaxonRegistry registry) {
        String newick(int index) { return trees.get(index).toNewick(registry); }
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
