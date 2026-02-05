import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import preprocessing.GeneTrees;
import utils.Config;
import core.InferenceDP;
import core.SpeciesTreeScorer;
import tree.RangeBipartition;
import tree.MixedBipartition;
import tree.Tree;

/**
 * Main entry point for phylogeny project with GeneTrees processing.
 * 
 * This implementation uses the actual GeneTrees class to parse and process
 * gene trees from Newick format input files.
 */
public class Main {

    /**
     * Main method that handles command line arguments and orchestrates the
     * analysis.
     */
    public static void main(String[] args) throws IOException {

        String inputFilePath = null;
        String outputFilePath = null;
        String computationMode = null;
        String scoreType = null;
        // expansionMethod and distanceMethod are fixed in this branch
        boolean verboseExpansion = false;
        String branchSupport = null;
        double lambda = 0.5;
        boolean useMixedBipartitions = false; // Cross-tree recombination flag (default: OFF)
        String speciesTreePath = null; // For score-only mode

        // Parse command line arguments
        for (int i = 0; i < args.length; i++) {
            if ((args[i].equals("-i") || args[i].equals("--input")) && i + 1 < args.length) {
                inputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            } else if ((args[i].equals("-o") || args[i].equals("--output")) && i + 1 < args.length) {
                outputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            } else if ((args[i].equals("-c") || args[i].equals("--score") || args[i].equals("--species-tree"))
                    && i + 1 < args.length) {
                speciesTreePath = args[i + 1];
                i++; // Skip next argument as it's the species tree path
            } else if (args[i].equals("--cpu")) {
                computationMode = "CPU_SINGLE";
            } else if (args[i].equals("--cpu-parallel")) {
                computationMode = "CPU_PARALLEL";
            } else if (args[i].equals("--gpu") || args[i].equals("--gpu-parallel")) {
                computationMode = "GPU_PARALLEL";
            } else if (args[i].equals("--triplet")) {
                scoreType = "TRIPLET";
            } else if (args[i].equals("--quartet")) {
                scoreType = "QUARTET";
            } else if (args[i].equals("--opt") && i + 1 < args.length) {
                scoreType = args[i + 1].toUpperCase();
                i++;
            } else if (args[i].equals("--score-type") && i + 1 < args.length) {
                scoreType = args[i + 1].toUpperCase();
                i++;
            } else if ((args[i].equals("-m") || args[i].equals("--mode")) && i + 1 < args.length) {
                computationMode = args[i + 1];
                i++; // Skip next argument as it's the mode
            } else if (args[i].equals("--expansion") || args[i].equals("-e")) {
                useMixedBipartitions = true;
                if (i + 1 < args.length && !args[i + 1].startsWith("-")) {
                    // legacy value (ignored)
                    i++;
                }
            } else if (args[i].equals("-v") || args[i].equals("--verbose")) {
                verboseExpansion = true;
            } else if ((args[i].equals("-s") || args[i].equals("--support") || args[i].equals("--branch-support"))
                    && i + 1 < args.length) {
                branchSupport = args[i + 1];
                i++; // Skip next argument as it's the support type
            } else if (args[i].equals("--lambda") && i + 1 < args.length) {
                try {
                    lambda = Double.parseDouble(args[i + 1]);
                    i++; // Skip next argument as it's the lambda value
                } catch (NumberFormatException e) {
                    System.err.println("Error: Invalid lambda value '" + args[i + 1] + "'");
                    System.exit(-1);
                }
            } else if (args[i].equals("--use-mixed") || args[i].equals("--extend-candidates")
                    || args[i].equals("--mixed") || args[i].equals("--no-mixed")) {
                System.err.println(
                        "Error: Mixed bipartitions are enabled only with --expansion. Remove deprecated mixed flags.");
                System.exit(-1);
            }
        }

        // Validate required arguments
        // Score-only mode: -i gene_trees.tre -c species_tree.tre (no -o needed)
        // Inference mode: -i gene_trees.tre -o output.tre
        if (inputFilePath == null || (outputFilePath == null && speciesTreePath == null)) {
            System.out.println("Usage:");
            System.out.println("  Inference mode: java Main -i <gene_trees> -o <output_file> [options]");
            System.out.println("  Score mode:     java Main -i <gene_trees> -c <species_tree> [options]");
            System.out.println("");
            System.out.println("Options:");
            System.out.println(
                    "  -c <tree>            Calculate triplet score between gene trees and given species tree");
            System.out.println("  --score <tree>       Same as -c");
            System.out.println("  --species-tree <tree> Same as -c");
            System.out.println("  --cpu                Computation mode: CPU_SINGLE");
            System.out.println("  --cpu-parallel        Computation mode: CPU_PARALLEL");
            System.out.println("  --gpu                Computation mode: GPU_PARALLEL");
            System.out.println("  -m, --mode <mode>     Computation mode: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
            System.out.println("  --score-type <type>   Score type: TRIPLET or QUARTET (default: TRIPLET)");
            System.out.println("  --triplet             Same as --score-type TRIPLET");
            System.out.println("  --quartet             Same as --score-type QUARTET");
            System.out.println("  --opt <type>          Alias for --score-type (TRIPLET or QUARTET)");
            System.out.println("  --expansion, -e        Enable mixed bipartitions (default: OFF)");
            System.out.println(
                    "  -s, --support <type>  Branch support: NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL");
            System.out.println("  --lambda <val>        Lambda parameter for branch support (default: 0.5)");
            System.out.println("  -v, --verbose         Verbose expansion output");
            System.out.println("  (Mixed bipartitions are enabled when --expansion is set)");
            System.exit(-1);
        }

        // Validate input file exists
        File inputFile = new File(inputFilePath);
        if (!inputFile.exists()) {
            System.err.println("Error: Input file '" + inputFilePath + "' does not exist.");
            System.exit(-1);
        }

        // Set computation mode if specified
        if (computationMode != null) {
            try {
                Config.COMPUTATION_MODE = Config.ComputationMode.valueOf(computationMode);
            } catch (IllegalArgumentException e) {
                System.err.println("Error: Invalid computation mode '" + computationMode + "'");
                System.err.println("Valid modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
                System.exit(-1);
            }
        }
        if (scoreType != null) {
            try {
                Config.SCORE_TYPE = Config.ScoreType.valueOf(scoreType);
            } catch (IllegalArgumentException e) {
                System.err.println("Error: Invalid score type '" + scoreType + "'");
                System.err.println("Valid types: TRIPLET, QUARTET");
                System.exit(-1);
            }
        }

        // Expansion method is fixed to NONE in this branch
        utils.BipartitionExpansionConfig.EXPANSION_METHOD = utils.BipartitionExpansionConfig.ExpansionMethod.NONE;

        // Set verbose expansion if specified
        if (verboseExpansion) {
            utils.BipartitionExpansionConfig.VERBOSE_EXPANSION = true;
            System.out.println("Verbose expansion output enabled.");
        }

        // Mixed bipartitions are enabled by --expansion

        // Determine mode
        boolean scoreMode = (speciesTreePath != null);

        System.out.println("Input file: " + inputFilePath);
        if (scoreMode) {
            System.out.println("Mode: SCORE (calculate triplet score for given species tree)");
            System.out.println("Species tree: " + speciesTreePath);
        } else {
            System.out.println("Mode: INFERENCE (find optimal species tree)");
            System.out.println("Output file: " + outputFilePath);
        }
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Score type: " + Config.SCORE_TYPE);
        if (!scoreMode) {
            System.out.println("Expansion method: " + utils.BipartitionExpansionConfig.EXPANSION_METHOD);
            if (utils.BipartitionExpansionConfig.isDistanceExpansionEnabled()) {
                System.out.println("Distance method: " + utils.BipartitionExpansionConfig.DISTANCE_METHOD);
            }
            System.out.println(
                    "Cross-tree recombination (via --expansion): " + (useMixedBipartitions ? "ENABLED" : "disabled"));
            if (branchSupport != null) {
                System.out.println("Branch support: " + branchSupport);
                System.out.println("Lambda parameter: " + lambda);
            }
        }

        long startTime = System.nanoTime();

        // Process gene trees
        GeneTrees geneTrees = new GeneTrees(inputFilePath);
        geneTrees.readTaxaNames(); // Ensure taxaMap is initialized
        geneTrees.readGeneTrees(null);

        // ================================================================
        // SCORE MODE: Calculate triplet score for given species tree
        // ================================================================
        if (scoreMode) {
            // Validate species tree file exists
            File speciesFile = new File(speciesTreePath);
            if (!speciesFile.exists()) {
                System.err.println("Error: Species tree file '" + speciesTreePath + "' does not exist.");
                System.exit(-1);
            }

            // Read species tree newick
            String speciesNewick = null;
            try (BufferedReader reader = new BufferedReader(new FileReader(speciesTreePath))) {
                speciesNewick = reader.readLine();
                if (speciesNewick != null) {
                    speciesNewick = speciesNewick.trim();
                }
            }

            if (speciesNewick == null || speciesNewick.isEmpty()) {
                System.err.println("Error: Species tree file is empty.");
                System.exit(-1);
            }

            System.out.println("\nParsing species tree...");
            Tree speciesTree = new Tree(speciesNewick, geneTrees.taxaMap);
            System.out.println("Species tree parsed: " + speciesTree.leavesCount + " leaves");

            // Calculate score (triplet or quartet)
            SpeciesTreeScorer scorer = new SpeciesTreeScorer(geneTrees);
            double score = scorer.calculateScore(speciesTree);

            long endTime = System.nanoTime();
            double duration = (endTime - startTime) / 1_000_000_000.0;

            System.out.println("\n========================================");
            if (Config.SCORE_TYPE == Config.ScoreType.QUARTET) {
                System.out.println("QUARTET_SCORE: " + score);
            } else {
                System.out.println("TRIPLET_SCORE: " + score);
            }

            // Calculate normalized score
            // Triplet: score / (k * (n choose 3))
            // Quartet: score / (k * (n choose 4))
            double k = (double) geneTrees.geneTrees.size();
            double n = (double) speciesTree.leavesCount;

            if (Config.SCORE_TYPE == Config.ScoreType.QUARTET) {
                if (n >= 4) {
                    double maxQuartetsPerTree = (n * (n - 1) * (n - 2) * (n - 3)) / 24.0;
                    double maxPossibleScore = k * maxQuartetsPerTree;
                    double normalizedScore = score / maxPossibleScore;
                    System.out.println("NORMALIZED_QUARTET_SCORE: " + normalizedScore);
                } else {
                    System.out.println("NORMALIZED_QUARTET_SCORE: Undefined (n < 4)");
                }
            } else {
                if (n >= 3) {
                    double maxTripletsPerTree = (n * (n - 1) * (n - 2)) / 6.0;
                    double maxPossibleScore = k * maxTripletsPerTree;
                    double normalizedScore = score / maxPossibleScore;
                    System.out.println("NORMALIZED_TRIPLET_SCORE: " + normalizedScore);
                } else {
                    System.out.println("NORMALIZED_TRIPLET_SCORE: Undefined (n < 3)");
                }
            }

            System.out.println("========================================");
            System.out.println("Time taken: " + duration + " seconds");
            System.out.println("Score calculation completed successfully!");

            return;
        }

        // ================================================================
        // INFERENCE MODE: Find optimal species tree via DP
        // ================================================================

        // Generate candidate bipartitions with cross-tree recombination extension
        System.out.println("Generating candidate bipartitions...");

        // Generate mixed bipartitions via cross-tree recombination
        // These will only be used in DP if expansion is enabled
        List<RangeBipartition> candidates = geneTrees.generateExtendedCandidateBipartitions(useMixedBipartitions);
        System.out.println("Total candidate bipartitions (gene tree): " + candidates.size());

        // Report mixed bipartitions generated by cross-tree recombination
        List<MixedBipartition> mixedBips = geneTrees.getMixedBipartitions();
        if (mixedBips != null && !mixedBips.isEmpty()) {
            System.out.println("Mixed bipartitions from cross-tree recombination: " + mixedBips.size());
            long crossTree = mixedBips.stream().filter(MixedBipartition::isCrossTree).count();
            System.out.println("  - Cross-tree (truly new): " + crossTree);
            System.out.println("  - Same-tree: " + (mixedBips.size() - crossTree));
        }

        // Run inference
        InferenceDP inference = new InferenceDP(geneTrees, candidates);

        // Enable mixed bipartitions in DP if expansion is enabled
        if (useMixedBipartitions && mixedBips != null && !mixedBips.isEmpty()) {
            System.out.println("\nEnabling mixed bipartitions in DP inference...");
            inference.enableMixedBipartitions(mixedBips);
        }

        double score = inference.solve();
        Tree resultTree = inference.reconstructTree();

        // Calculate branch support if requested
        if (branchSupport != null && !branchSupport.equals("NONE")) {
            System.out.println("\nCalculating branch support...");

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

            core.BranchSupportCalculator supportCalculator = new core.BranchSupportCalculator(geneTrees, resultTree,
                    lambda, annotationType);

            // Validate quartet frequencies for debugging (optional)
            if (verboseExpansion) {
                supportCalculator.validateQuartetFrequencies();
            }

            // Annotate branches
            supportCalculator.annotateBranches();

            // Print statistics
            core.BranchSupportCalculator.BranchSupportStatistics stats = supportCalculator.calculateStatistics();
            System.out.println("\n" + stats.toString());
        }

        // Write output
        try (FileWriter writer = new FileWriter(outputFilePath)) {
            writer.write(resultTree.getNewickFormat());
        }

        long endTime = System.nanoTime();
        double duration = (endTime - startTime) / 1_000_000_000.0; // Convert to seconds

        System.out.println("\n========================================");
        if (Config.SCORE_TYPE == Config.ScoreType.QUARTET) {
            System.out.println("OPTIMAL_QUARTET_SCORE: " + score);
        } else {
            System.out.println("OPTIMAL_TRIPLET_SCORE: " + score);
        }
        System.out.println("========================================");
        System.out.println("Time taken: " + duration + " seconds");
        System.out.println("Program completed successfully!");
        System.out.println("Output written to: " + outputFilePath);
    }

    /**
     * Processes gene trees using the GeneTrees class and returns analysis results.
     * 
     * @param inputFilePath Path to the input file containing gene trees in Newick
     *                      format
     * @return Formatted string with analysis results
     * @throws FileNotFoundException if the input file cannot be read
     */
    private static String processGeneTrees(String inputFilePath) throws FileNotFoundException {
        System.out.println("Initializing GeneTrees...");

        // Create GeneTrees object and read taxa names
        GeneTrees geneTrees = new GeneTrees(inputFilePath);
        var taxaMap = geneTrees.readTaxaNames();

        System.out.println("Reading and parsing gene trees...");

        // Read and process all gene trees
        geneTrees.readGeneTrees(null); // No distance matrix needed for basic analysis

        // Debug output
        // debugOutput(geneTrees);

        System.out.println(geneTrees.geneTrees.get(0).isRooted);

        // Test InferenceDP algorithm
        System.out.println("Testing InferenceDP algorithm...");
        List<RangeBipartition> candidates = new ArrayList<>(geneTrees.rangeBipartitions.keySet());

        if (!candidates.isEmpty()) {
            InferenceDP dp = new InferenceDP(geneTrees, candidates);
            double maxScore = dp.solve();

            System.out.println("DP Algorithm completed with maximum score: " + maxScore);

            Tree reconstructedTree = dp.reconstructTree();
            if (reconstructedTree != null && reconstructedTree.root != null) {
                System.out.println("Tree reconstruction successful");
                return reconstructedTree.getNewickFormat();
            }
        }

        return "";
    }

    /**
     * Debug function that prints detailed analysis information to console
     */
    private static void debugOutput(GeneTrees geneTrees) {
        // Gather analysis information
        int geneTreeCount = geneTrees.geneTrees.size();
        int taxaCount = geneTrees.realTaxaCount;
        // int uniquePartitions = geneTrees.triPartitions.size();
        int uniqueRangeBipartitions = geneTrees.rangeBipartitions.size();

        System.out.println("Processing complete:");
        System.out.println("  - Gene trees processed: " + geneTreeCount);
        System.out.println("  - Taxa found: " + taxaCount);
        // System.out.println(" - Unique tripartitions: " + uniquePartitions);
        System.out.println("  - Unique RangeBipartitions: " + uniqueRangeBipartitions);

        // Print taxa names
        System.out.print("Taxa names: ");
        for (int i = 0; i < geneTrees.taxonIdToLabel.length; i++) {
            if (i > 0)
                System.out.print(", ");
            System.out.print(geneTrees.taxonIdToLabel[i]);
        }
        System.out.println();

        // Print RangeBipartitions with counts
        System.out.println("RangeBipartitions:");
        for (var entry : geneTrees.rangeBipartitions.entrySet()) {
            System.out.println("  " + entry.getKey().toString() + " : " + entry.getValue());
        }
    }

    /**
     * Writes analysis results to the specified output file.
     * 
     * @param outputFilePath Path to the output file
     * @param content        Content to write to the file
     * @throws IOException if there's an error writing to the file
     */
    private static void writeResults(String outputFilePath, String content) throws IOException {
        FileWriter writer = new FileWriter(outputFilePath);
        writer.write(content);
        writer.close();
    }
}
