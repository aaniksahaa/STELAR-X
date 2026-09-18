package stelarx;

/** Focused, dependency-free tests for the user-facing CLI mappings. */
public final class CliPresetsTest {
    private record SearchExpected(
        Config.SearchMode mode,
        boolean autocomplete,
        boolean consensus,
        boolean quadratic,
        boolean randomResolution,
        boolean largePolytomies,
        boolean inputPolytomies) {}

    public static void main(String[] args) {
        Config cfg = Config.getInstance();
        check(!cfg.isKeepPolytomyDuringInference(),
            "default inference polytomy refinement");
        cfg.setKeepPolytomyDuringInference(true);
        check(cfg.isKeepPolytomyDuringInference(),
            "keep-polytomy-during-inference config");
        cfg.setKeepPolytomyDuringInference(false);
        SearchExpected[] expected = {
            new SearchExpected(Config.SearchMode.LOCAL, false, false, false, false, false, false),
            new SearchExpected(Config.SearchMode.FULL,  true,  false, false, false, false, false),
            new SearchExpected(Config.SearchMode.FULL,  true,  true,  true,  true,  true,  true)
        };

        // Descending order proves that every lower preset resets higher-level flags.
        for (int level = 3; level >= 1; level--) {
            cfg.setCompletionMethod(Config.CompletionMethod.DISTANCE);
            cfg.setStepBFastRestriction(false);
            CliPresets.applySearchSpace("S" + level, cfg);
            assertSearch("S" + level, cfg, expected[level - 1]);
            check(cfg.getCompletionMethod() == Config.CompletionMethod.SIMILARITY,
                "S" + level + " completion method");
            check(cfg.isStepBFastRestriction(), "S" + level + " fast restriction");
            CliPresets.applySearchSpace(Integer.toString(level), cfg);
            assertSearch(Integer.toString(level), cfg, expected[level - 1]);
        }

        Config.WeightIntersectionMethod[] methods = {
            Config.WeightIntersectionMethod.SMALLER_SIDE_TRAVERSAL,
            Config.WeightIntersectionMethod.PREFIX_SUM,
            Config.WeightIntersectionMethod.SIMPLE_TREE_WALK,
            Config.WeightIntersectionMethod.BITSET
        };
        for (int level = 1; level <= 4; level++) {
            CliPresets.applyIntersectionMethod("I" + level, cfg);
            check(cfg.getWeightIntersectionMethod() == methods[level - 1], "I" + level);
            CliPresets.applyIntersectionMethod(Integer.toString(level), cfg);
            check(cfg.getWeightIntersectionMethod() == methods[level - 1], Integer.toString(level));
        }

        // Legacy names remain accepted through the same mapping.
        CliPresets.applyIntersectionMethod("smaller-side-traversal", cfg);
        check(cfg.getWeightIntersectionMethod() == methods[0], "legacy smaller-side name");
        CliPresets.applyIntersectionMethod("prefix_sum", cfg);
        check(cfg.getWeightIntersectionMethod() == methods[1], "legacy prefix name");
        CliPresets.applyIntersectionMethod("simple-tree-walk", cfg);
        check(cfg.getWeightIntersectionMethod() == methods[2], "legacy tree-walk name");
        CliPresets.applyIntersectionMethod("bitset", cfg);
        check(cfg.getWeightIntersectionMethod() == methods[3], "legacy bitset name");

        expectInvalid(() -> CliPresets.applySearchSpace("S4", cfg));
        expectInvalid(() -> CliPresets.applySearchSpace("S8", cfg));
        expectInvalid(() -> CliPresets.applyIntersectionMethod("I5", cfg));
        System.out.println("CLI preset mappings: PASS");
    }

    private static void assertSearch(String label, Config cfg, SearchExpected e) {
        check(cfg.getSearchMode() == e.mode, label + " search mode");
        check(cfg.isAutoCompleteIncompleteTrees() == e.autocomplete, label + " autocomplete");
        check(cfg.isConsensusExperimental() == e.consensus, label + " consensus");
        check(cfg.isStepBQuadraticNnBalls() == e.quadratic, label + " quadratic");
        check(cfg.isStepBRandomLeftoverResolution() == e.randomResolution, label + " random resolution");
        check(cfg.isStepBProcessLargePolytomies() == e.largePolytomies, label + " large polytomies");
        check(cfg.isResolveInputGeneTreePolytomies() == e.inputPolytomies, label + " input polytomies");
    }

    private static void expectInvalid(Runnable action) {
        try {
            action.run();
            throw new AssertionError("expected IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
            // Expected.
        }
    }

    private static void check(boolean condition, String label) {
        if (!condition) throw new AssertionError("failed: " + label);
    }
}
