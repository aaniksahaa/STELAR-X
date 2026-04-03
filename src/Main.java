import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;

import preprocessing.GeneTrees;
import utils.Config;
import utils.PhaseLogger;
import core.SpeciesTreeScorer;
import tree.Tree;
import taxon.Taxon;
import taxon.TaxonRegistry;
import hash.TaxonHasher;
import hash.PrefixHashArrays;
import cluster.ClusterTable;
import partition.PartitionTable;
import dp.DPTable;
import dp.Inference;
import weight.WeightTable;

/**
 * Main entry point for STELAR-X phylogeny inference.
 *
 * Uses the ASTRAL-X-style pipeline to dramatically reduce peak RAM:
 *   Phase 1  Parse    -> trees, registry
 *   Phase 2  Hash     -> hasher, pref         (hasher=null after)
 *   Phase 3  Clusters -> clusterTable         (uses pref)
 *   Phase 4  Partitions -> partTable          (uses pref)
 *   Phase 5  DP space   -> dpTable            (uses pref)
 *            *** pref=null; System.gc(); ***  <- THE KEY MEMORY WIN
 *   Phase 6  Weights  -> weightTable          (uses Tree.postorderArray/positionMap)
 *   Phase 7  Inference -> species tree        (uses dpTable + weightTable)
 */
public class Main {

    public static void main(String[] args) throws IOException {

        String inputFilePath   = null;
        String outputFilePath  = null;
        String computationMode = null;
        boolean verboseExpansion = false;
        String branchSupport   = null;
        double lambda          = 0.5;
        boolean useMixedBipartitions = false;
        String speciesTreePath = null;

        // ── Parse CLI arguments ────────────────────────────────────────────────
        for (int i = 0; i < args.length; i++) {
            if ((args[i].equals("-i") || args[i].equals("--input")) && i + 1 < args.length) {
                inputFilePath = args[++i];
            } else if ((args[i].equals("-o") || args[i].equals("--output")) && i + 1 < args.length) {
                outputFilePath = args[++i];
            } else if ((args[i].equals("-c") || args[i].equals("--score") || args[i].equals("--species-tree"))
                    && i + 1 < args.length) {
                speciesTreePath = args[++i];
            } else if (args[i].equals("--cpu")) {
                computationMode = "CPU_SINGLE";
            } else if (args[i].equals("--cpu-parallel")) {
                computationMode = "CPU_PARALLEL";
            } else if (args[i].equals("--gpu") || args[i].equals("--gpu-parallel")) {
                computationMode = "GPU_PARALLEL";
            } else if ((args[i].equals("-m") || args[i].equals("--mode")) && i + 1 < args.length) {
                computationMode = args[++i];
            } else if (args[i].equals("--expansion") || args[i].equals("-e")) {
                useMixedBipartitions = true;
                if (i + 1 < args.length && !args[i + 1].startsWith("-")) {
                    i++; // legacy value – ignored
                }
            } else if (args[i].equals("-v") || args[i].equals("--verbose")) {
                verboseExpansion = true;
            } else if ((args[i].equals("-s") || args[i].equals("--support") || args[i].equals("--branch-support"))
                    && i + 1 < args.length) {
                branchSupport = args[++i];
            } else if (args[i].equals("--lambda") && i + 1 < args.length) {
                try {
                    lambda = Double.parseDouble(args[++i]);
                } catch (NumberFormatException e) {
                    System.err.println("Error: Invalid lambda value '" + args[i] + "'");
                    System.exit(-1);
                }
            } else if (args[i].equals("--use-mixed") || args[i].equals("--extend-candidates")
                    || args[i].equals("--mixed") || args[i].equals("--no-mixed")) {
                System.err.println(
                        "Error: Mixed bipartitions are enabled only with --expansion. Remove deprecated mixed flags.");
                System.exit(-1);
            }
        }

        // ── Validate arguments ─────────────────────────────────────────────────
        if (inputFilePath == null || (outputFilePath == null && speciesTreePath == null)) {
            System.out.println("Usage:");
            System.out.println("  Inference mode: java Main -i <gene_trees> -o <output_file> [options]");
            System.out.println("  Score mode:     java Main -i <gene_trees> -c <species_tree> [options]");
            System.out.println("");
            System.out.println("Options:");
            System.out.println("  -c <tree>             Calculate triplet score between gene trees and given species tree");
            System.out.println("  --score <tree>        Same as -c");
            System.out.println("  --species-tree <tree> Same as -c");
            System.out.println("  --cpu                 Computation mode: CPU_SINGLE");
            System.out.println("  --cpu-parallel        Computation mode: CPU_PARALLEL");
            System.out.println("  --gpu                 Computation mode: GPU_PARALLEL");
            System.out.println("  -m, --mode <mode>     Computation mode: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
            System.out.println("  --expansion, -e       Enable mixed bipartitions (default: OFF)");
            System.out.println("  -s, --support <type>  Branch support: NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL");
            System.out.println("  --lambda <val>        Lambda parameter for branch support (default: 0.5)");
            System.out.println("  -v, --verbose         Verbose expansion output");
            System.exit(-1);
        }

        File inputFile = new File(inputFilePath);
        if (!inputFile.exists()) {
            System.err.println("Error: Input file '" + inputFilePath + "' does not exist.");
            System.exit(-1);
        }

        if (computationMode != null) {
            try {
                Config.COMPUTATION_MODE = Config.ComputationMode.valueOf(computationMode);
            } catch (IllegalArgumentException e) {
                System.err.println("Error: Invalid computation mode '" + computationMode + "'");
                System.err.println("Valid modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
                System.exit(-1);
            }
        }

        utils.BipartitionExpansionConfig.EXPANSION_METHOD =
                utils.BipartitionExpansionConfig.ExpansionMethod.NONE;

        if (verboseExpansion) {
            utils.BipartitionExpansionConfig.VERBOSE_EXPANSION = true;
            System.out.println("Verbose expansion output enabled.");
        }

        boolean scoreMode = (speciesTreePath != null);

        printRuntimeEnvironment(Config.COMPUTATION_MODE);

        System.out.println("Input file: " + inputFilePath);
        if (scoreMode) {
            System.out.println("Mode: SCORE (calculate triplet score for given species tree)");
            System.out.println("Species tree: " + speciesTreePath);
        } else {
            System.out.println("Mode: INFERENCE (find optimal species tree)");
            System.out.println("Output file: " + outputFilePath);
            System.out.println("Expansion method: " + utils.BipartitionExpansionConfig.EXPANSION_METHOD);
            System.out.println("Cross-tree recombination (via --expansion): " +
                    (useMixedBipartitions ? "ENABLED" : "disabled"));
            if (branchSupport != null) {
                System.out.println("Branch support: " + branchSupport);
                System.out.println("Lambda parameter: " + lambda);
            }
        }

        long startTime = System.nanoTime();

        // ================================================================
        // SCORE MODE: delegate to SpeciesTreeScorer (unchanged path)
        // ================================================================
        if (scoreMode) {
            File speciesFile = new File(speciesTreePath);
            if (!speciesFile.exists()) {
                System.err.println("Error: Species tree file '" + speciesTreePath + "' does not exist.");
                System.exit(-1);
            }

            String speciesNewick = null;
            try (BufferedReader reader = new BufferedReader(new FileReader(speciesTreePath))) {
                speciesNewick = reader.readLine();
                if (speciesNewick != null) speciesNewick = speciesNewick.trim();
            }
            if (speciesNewick == null || speciesNewick.isEmpty()) {
                System.err.println("Error: Species tree file is empty.");
                System.exit(-1);
            }

            GeneTrees geneTrees = new GeneTrees(inputFilePath);
            geneTrees.readTaxaNames();
            geneTrees.readGeneTrees(null);

            System.out.println("\nParsing species tree...");
            Tree speciesTree = new Tree(speciesNewick, geneTrees.taxaMap);
            System.out.println("Species tree parsed: " + speciesTree.leavesCount + " leaves");

            SpeciesTreeScorer scorer = new SpeciesTreeScorer(geneTrees);
            double score = scorer.calculateScore(speciesTree);

            long endTime = System.nanoTime();
            double duration = (endTime - startTime) / 1_000_000_000.0;

            System.out.println("\n========================================");
            System.out.println("TRIPLET_SCORE: " + score);

            double k = (double) geneTrees.geneTrees.size();
            double n = (double) speciesTree.leavesCount;
            if (n >= 3) {
                double maxPossibleScore = k * (n * (n - 1) * (n - 2)) / 6.0;
                System.out.println("NORMALIZED_TRIPLET_SCORE: " + score / maxPossibleScore);
            } else {
                System.out.println("NORMALIZED_TRIPLET_SCORE: Undefined (n < 3)");
            }
            System.out.println("========================================");
            System.out.println("Time taken: " + duration + " seconds");
            System.out.println("Score calculation completed successfully!");
            return;
        }

        // ================================================================
        // INFERENCE MODE: ASTRAL-X-style 7-phase pipeline
        // ================================================================

        // ── Phase 1: Parse ─────────────────────────────────────────────────────
        long t1 = PhaseLogger.begin("Phase 1  Parse", false);

        GeneTrees geneTrees = new GeneTrees(inputFilePath);
        geneTrees.readTaxaNames();
        geneTrees.readGeneTrees(null);

        List<tree.Tree> trees = geneTrees.geneTrees;
        int numTaxa = geneTrees.realTaxaCount;

        // Assign treeIndex to each tree (needed by pipeline classes)
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).treeIndex = i;
        }

        // Build TaxonRegistry from geneTrees.taxaMap  (name -> Taxon)
        TaxonRegistry registry = buildRegistry(geneTrees.taxaMap, geneTrees.taxonIdToLabel, numTaxa);

        PhaseLogger.end("Phase 1  Parse", t1, false);

        // ── Phase 2: Hash ──────────────────────────────────────────────────────
        long t2 = PhaseLogger.begin("Phase 2  Hash", false);

        TaxonHasher hasher = new TaxonHasher(numTaxa, 2, 0xDEADBEEFL);
        PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);
        hasher = null; // free after pref is built

        PhaseLogger.end("Phase 2  Hash", t2, false);

        // ── Phase 3: Clusters ──────────────────────────────────────────────────
        long t3 = PhaseLogger.begin("Phase 3  Clusters", false);

        ClusterTable clusterTable = new ClusterTable(trees, pref, numTaxa);

        PhaseLogger.end("Phase 3  Clusters", t3, false);

        // ── Phase 4: Partitions ────────────────────────────────────────────────
        long t4 = PhaseLogger.begin("Phase 4  Partitions", false);

        PartitionTable partTable = new PartitionTable(trees, pref);

        PhaseLogger.end("Phase 4  Partitions", t4, false);

        // ── Phase 5: DP space ──────────────────────────────────────────────────
        long t5 = PhaseLogger.begin("Phase 5  DP space", false);

        DPTable dpTable = new DPTable(trees, pref, clusterTable);

        PhaseLogger.end("Phase 5  DP space", t5, false);

        // *** THE KEY MEMORY WIN: free prefix arrays before weight calc ***
        pref = null;
        long gcBefore = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        System.gc();
        long gcAfter = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        System.out.println("Pre-Phase-6 GC: heap " + gcBefore / 1_000_000 + " MB -> "
                + gcAfter / 1_000_000 + " MB (freed " + (gcBefore - gcAfter) / 1_000_000 + " MB)");

        // ── Phase 6: Weights ───────────────────────────────────────────────────
        long t6 = PhaseLogger.begin("Phase 6  Weights", false);

        WeightTable weightTable = new WeightTable(dpTable, partTable, clusterTable, trees);

        PhaseLogger.end("Phase 6  Weights", t6, false);

        // ── Phase 7: Inference ─────────────────────────────────────────────────
        long t7 = PhaseLogger.begin("Phase 7  Inference", false);

        Inference inference = new Inference();
        String newick = inference.run(dpTable, weightTable, clusterTable, trees, registry);

        PhaseLogger.end("Phase 7  Inference", t7, false);

        // ── Write output ───────────────────────────────────────────────────────

        // Branch support annotation (if requested) requires a Tree object.
        // Parse the Newick string back into a Tree for annotation.
        if (branchSupport != null && !branchSupport.equals("NONE")) {
            System.out.println("\nCalculating branch support...");

            Tree resultTree = new Tree(newick, geneTrees.taxaMap);

            core.BranchSupportCalculator.BranchAnnotationType annotationType;
            try {
                switch (branchSupport.toUpperCase()) {
                    case "POSTERIOR":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.POSTERIOR_ONLY;
                        break;
                    case "DETAILED":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.DETAILED;
                        break;
                    case "LENGTH":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.BRANCH_LENGTH_ONLY;
                        break;
                    case "BOTH":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.POSTERIOR_AND_LENGTH;
                        break;
                    case "PVALUE":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.PVALUE_ONLY;
                        break;
                    case "ALL":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.ALL;
                        break;
                    default:
                        System.err.println("Error: Invalid branch support type '" + branchSupport + "'");
                        System.err.println("Valid types: NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL");
                        System.exit(-1);
                        return;
                }
            } catch (Exception e) {
                System.err.println("Error parsing branch support type: " + e.getMessage());
                System.exit(-1);
                return;
            }

            core.BranchSupportCalculator supportCalculator =
                    new core.BranchSupportCalculator(geneTrees, resultTree, lambda, annotationType);

            if (verboseExpansion) {
                supportCalculator.validateQuartetFrequencies();
            }
            supportCalculator.annotateBranches();

            core.BranchSupportCalculator.BranchSupportStatistics stats =
                    supportCalculator.calculateStatistics();
            System.out.println("\n" + stats.toString());

            newick = resultTree.getNewickFormat();
        }

        try (FileWriter writer = new FileWriter(outputFilePath)) {
            writer.write(newick);
        }

        long endTime = System.nanoTime();
        double duration = (endTime - startTime) / 1_000_000_000.0;

        System.out.println("\n========================================");
        System.out.println("OPTIMAL_TRIPLET_SCORE: " + inference.getDPScore(dpTable.getRootHash()));
        System.out.println("========================================");
        System.out.println("Time taken: " + duration + " seconds");
        System.out.println("Program completed successfully!");
        System.out.println("Output written to: " + outputFilePath);
    }

    // ── Helpers ────────────────────────────────────────────────────────────────

    /**
     * Build a TaxonRegistry from the taxaMap produced by GeneTrees.
     * Uses the taxonIdToLabel array to ensure IDs match exactly.
     */
    private static TaxonRegistry buildRegistry(Map<String, Taxon> taxaMap,
                                               String[] taxonIdToLabel,
                                               int numTaxa) {
        TaxonRegistry registry = new TaxonRegistry();
        // Register names in ID order so getId(name) == taxon.id
        if (taxonIdToLabel != null) {
            for (int id = 0; id < numTaxa && id < taxonIdToLabel.length; id++) {
                String name = taxonIdToLabel[id];
                if (name != null) {
                    registry.register(name);
                }
            }
        } else {
            // Fallback: register in whatever map order
            for (String name : taxaMap.keySet()) {
                registry.register(name);
            }
        }
        registry.lock();
        return registry;
    }

    /** Prints a concise runtime environment summary. */
    private static void printRuntimeEnvironment(Config.ComputationMode mode) {
        int cpuCores      = Runtime.getRuntime().availableProcessors();
        long maxHeapBytes = Runtime.getRuntime().maxMemory();
        String heapStr    = formatBytes(maxHeapBytes);

        String deviceLine;
        switch (mode) {
            case GPU_PARALLEL:
                String gpuName = detectGpuName();
                deviceLine = (gpuName != null) ? "GPU \u2014 " + gpuName
                        : "GPU (device details unavailable \u2014 will attempt CUDA at runtime)";
                break;
            case CPU_PARALLEL:
                deviceLine = "CPU (" + cpuCores + " cores, multi-threaded)";
                break;
            case CPU_SINGLE:
            default:
                deviceLine = "CPU (single-threaded)";
                break;
        }

        String CYAN  = "\u001B[36m";
        String BOLD  = "\u001B[1m";
        String GREEN = "\u001B[32m";
        String RESET = "\u001B[0m";

        System.out.println();
        System.out.println(CYAN + "\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550" + RESET);
        System.out.println(BOLD + CYAN + "  STELAR-X  Runtime Environment" + RESET);
        System.out.println(CYAN + "\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550" + RESET);
        System.out.println("  Compute:  " + BOLD + GREEN + mode + RESET);
        System.out.println("  Device:   " + deviceLine);
        System.out.println("  Threads:  " + cpuCores + " available");
        System.out.println("  JVM Heap: " + heapStr);
        System.out.println(CYAN + "\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550" + RESET);
        System.out.println();
    }

    private static String detectGpuName() {
        try {
            Process proc = new ProcessBuilder(
                    "nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader,nounits")
                    .redirectErrorStream(true).start();
            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(proc.getInputStream()))) {
                String line = reader.readLine();
                int exitCode = proc.waitFor();
                if (exitCode == 0 && line != null && !line.trim().isEmpty()) {
                    String[] parts = line.split(",");
                    if (parts.length >= 2) {
                        String name   = parts[0].trim();
                        String memMB  = parts[1].trim();
                        try {
                            int mb = Integer.parseInt(memMB);
                            String memStr = (mb >= 1024)
                                    ? String.format("%.0f GB", mb / 1024.0) : (mb + " MB");
                            return name + " (" + memStr + ")";
                        } catch (NumberFormatException e) {
                            return name;
                        }
                    }
                    return line.trim();
                }
            }
        } catch (Exception ignored) {}
        return null;
    }

    private static String formatBytes(long bytes) {
        if (bytes >= 1024L * 1024 * 1024)
            return String.format("%.1f GB", bytes / (1024.0 * 1024 * 1024));
        else if (bytes >= 1024L * 1024)
            return String.format("%.0f MB", bytes / (1024.0 * 1024));
        else
            return bytes + " bytes";
    }
}
