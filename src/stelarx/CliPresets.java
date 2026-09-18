package stelarx;

import java.util.Locale;

/**
 * User-facing CLI presets.  These methods only translate friendly names into
 * the existing Config switches; they do not change inference algorithms.
 */
final class CliPresets {
    private CliPresets() {}

    static void applySearchSpace(String value, Config cfg) {
        String preset = normalize(value);

        // A preset is a complete, deterministic configuration.  This also makes
        // command-line ordering unsurprising: options are applied left-to-right,
        // so a later legacy flag can fine-tune a preset and a later preset resets
        // all search-space flags.
        cfg.setSearchMode(Config.SearchMode.LOCAL);
        cfg.setAutoCompleteIncompleteTrees(false);
        cfg.setCompletionMethod(Config.CompletionMethod.SIMILARITY);
        cfg.setConsensusExperimental(false);
        cfg.setStepBFastRestriction(true);
        cfg.setStepBQuadraticNnBalls(false);
        cfg.setStepBRandomLeftoverResolution(false);
        cfg.setStepBProcessLargePolytomies(false);
        cfg.setResolveInputGeneTreePolytomies(false);

        switch (preset) {
            case "1", "s1", "incomplete-local" -> {
                // Baseline: original trees and tree-local candidate splits.
            }
            case "2", "s2", "complete-full" ->
                applyCompleteFull(cfg);
            case "3", "s3", "exhaustive" -> {
                applyCompleteFull(cfg);
                cfg.setConsensusExperimental(true);
                cfg.setStepBQuadraticNnBalls(true);
                cfg.setStepBRandomLeftoverResolution(true);
                cfg.setStepBProcessLargePolytomies(true);
                cfg.setResolveInputGeneTreePolytomies(true);
            }
            default -> throw new IllegalArgumentException(
                "unknown search space '" + value + "' (expected S1-S3 or 1-3)");
        }
    }

    static void applyIntersectionMethod(String value, Config cfg) {
        String method = normalize(value);
        switch (method) {
            case "1", "i1", "smaller-side-traversal", "smaller-side", "smallerside", "legacy" ->
                cfg.setWeightIntersectionMethod(Config.WeightIntersectionMethod.SMALLER_SIDE_TRAVERSAL);
            case "2", "i2", "prefix-sum", "prefixsum", "prefix" ->
                cfg.setWeightIntersectionMethod(Config.WeightIntersectionMethod.PREFIX_SUM);
            case "3", "i3", "simple-tree-walk", "tree-walk", "treewalk", "simple" ->
                cfg.setWeightIntersectionMethod(Config.WeightIntersectionMethod.SIMPLE_TREE_WALK);
            case "4", "i4", "bitset", "bitsets", "bit-set" ->
                cfg.setWeightIntersectionMethod(Config.WeightIntersectionMethod.BITSET);
            default -> throw new IllegalArgumentException(
                "unknown intersection method '" + value + "' (expected I1-I4 or 1-4)");
        }
    }

    private static void applyCompleteFull(Config cfg) {
        cfg.setAutoCompleteIncompleteTrees(true);
        cfg.setSearchMode(Config.SearchMode.FULL);
    }

    private static String normalize(String value) {
        return value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
    }
}
