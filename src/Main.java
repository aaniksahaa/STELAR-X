import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import preprocessing.GeneTrees;
import utils.Config;
import core.InferenceDP;
import tree.STBipartition;
import tree.Tree;

/**
 * Main entry point for phylogeny project with GeneTrees processing.
 * 
 * This implementation uses the actual GeneTrees class to parse and process
 * gene trees from Newick format input files.
 */
public class Main {

    /**
     * Main method that handles command line arguments and orchestrates the analysis.
     */
    public static void main(String[] args) throws IOException {

        String inputFilePath = null;
        String outputFilePath = null;
        String computationMode = null;
        String expansionMethod = null;
        String distanceMethod = null;
        boolean verboseExpansion = false;
        boolean disableExpansion = false;
        String branchSupport = null;
        double lambda = 0.5;

        // Parse command line arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i") && i + 1 < args.length) {
                inputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            } else if (args[i].equals("-o") && i + 1 < args.length) {
                outputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            } else if (args[i].equals("-m") && i + 1 < args.length) {
                computationMode = args[i + 1];
                i++; // Skip next argument as it's the mode
            } else if (args[i].equals("-e") && i + 1 < args.length) {
                expansionMethod = args[i + 1];
                i++; // Skip next argument as it's the expansion method
            } else if (args[i].equals("-d") && i + 1 < args.length) {
                distanceMethod = args[i + 1];
                i++; // Skip next argument as it's the distance method
            } else if (args[i].equals("-v")) {
                verboseExpansion = true;
            } else if (args[i].equals("--no-expansion")) {
                disableExpansion = true;
            } else if (args[i].equals("-s") && i + 1 < args.length) {
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
            }
        }

        // Validate required arguments
        if (inputFilePath == null || outputFilePath == null) {
            System.out.println("Usage: java Main -i <input_file> -o <output_file> [options]");
            System.out.println("Options:");
            System.out.println("  -m <mode>     Computation mode: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
            System.out.println("  -e <method>   Expansion method: NONE, DISTANCE_ONLY, CONSENSUS_ONLY, DISTANCE_CONSENSUS, FULL");
            System.out.println("  -d <method>   Distance method: UPGMA, NEIGHBOR_JOINING, BOTH");
            System.out.println("  -s <support>  Branch support: NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL");
            System.out.println("  --lambda <val> Lambda parameter for branch support (default: 0.5)");
            System.out.println("  -v            Verbose expansion output");
            System.out.println("  --no-expansion Disable bipartition expansion");
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
        
        // Configure bipartition expansion
        if (disableExpansion) {
            utils.BipartitionExpansionConfig.EXPANSION_METHOD = utils.BipartitionExpansionConfig.ExpansionMethod.NONE;
        } else if (expansionMethod != null) {
            try {
                utils.BipartitionExpansionConfig.EXPANSION_METHOD = 
                    utils.BipartitionExpansionConfig.ExpansionMethod.valueOf(expansionMethod);
            } catch (IllegalArgumentException e) {
                System.err.println("Error: Invalid expansion method '" + expansionMethod + "'");
                System.err.println("Valid methods: NONE, DISTANCE_ONLY, CONSENSUS_ONLY, DISTANCE_CONSENSUS, FULL");
                System.exit(-1);
            }
        }
        
        // Set distance method if specified
        if (distanceMethod != null) {
            try {
                utils.BipartitionExpansionConfig.DISTANCE_METHOD = 
                    utils.BipartitionExpansionConfig.DistanceMethod.valueOf(distanceMethod);
            } catch (IllegalArgumentException e) {
                System.err.println("Error: Invalid distance method '" + distanceMethod + "'");
                System.err.println("Valid methods: UPGMA, NEIGHBOR_JOINING, BOTH");
                System.exit(-1);
            }
        }
        
        // Set verbose expansion if specified
        if (verboseExpansion) {
            utils.BipartitionExpansionConfig.VERBOSE_EXPANSION = true;
            System.out.println("Verbose expansion output enabled.");
        }

        System.out.println("Input file: " + inputFilePath);
        System.out.println("Output file: " + outputFilePath);
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Expansion method: " + utils.BipartitionExpansionConfig.EXPANSION_METHOD);
        if (utils.BipartitionExpansionConfig.isDistanceExpansionEnabled()) {
            System.out.println("Distance method: " + utils.BipartitionExpansionConfig.DISTANCE_METHOD);
        }
        if (branchSupport != null) {
            System.out.println("Branch support: " + branchSupport);
            System.out.println("Lambda parameter: " + lambda);
        }

        long startTime = System.nanoTime();

        // Process gene trees
        GeneTrees geneTrees = new GeneTrees(inputFilePath);
        geneTrees.readTaxaNames(); // Ensure taxaMap is initialized
        geneTrees.readGeneTrees(null);

        // Generate candidate bipartitions
        System.out.println("Generating candidate bipartitions...");
        List<STBipartition> candidates = geneTrees.generateCandidateBipartitions();
        System.out.println("Total candidate bipartitions: " + candidates.size());

        // Run inference
        InferenceDP inference = new InferenceDP(geneTrees, candidates);
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
            
            core.BranchSupportCalculator supportCalculator = 
                new core.BranchSupportCalculator(geneTrees, resultTree, lambda, annotationType);
            
            // Validate quartet frequencies for debugging (optional)
            if (verboseExpansion) {
                supportCalculator.validateQuartetFrequencies();
            }
            
            // Annotate branches
            supportCalculator.annotateBranches();
            
            // Print statistics
            core.BranchSupportCalculator.BranchSupportStatistics stats = 
                supportCalculator.calculateStatistics();
            System.out.println("\n" + stats.toString());
        }

        // Write output
        try (FileWriter writer = new FileWriter(outputFilePath)) {
            writer.write(resultTree.getNewickFormat());
        }

        long endTime = System.nanoTime();
        double duration = (endTime - startTime) / 1_000_000_000.0; // Convert to seconds

        System.out.println("Score: " + score);
        System.out.println("Time taken: " + duration + " seconds");
        System.out.println("Program completed successfully!");
        System.out.println("Output written to: " + outputFilePath);
    }

    /**
     * Processes gene trees using the GeneTrees class and returns analysis results.
     * 
     * @param inputFilePath Path to the input file containing gene trees in Newick format
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
        List<STBipartition> candidates = new ArrayList<>(geneTrees.stBipartitions.keySet());
        
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
        int uniqueSTBipartitions = geneTrees.stBipartitions.size();
        
        System.out.println("Processing complete:");
        System.out.println("  - Gene trees processed: " + geneTreeCount);
        System.out.println("  - Taxa found: " + taxaCount);
        // System.out.println("  - Unique tripartitions: " + uniquePartitions);
        System.out.println("  - Unique STBipartitions: " + uniqueSTBipartitions);
        
        // Print taxa names
        System.out.print("Taxa names: ");
        for (int i = 0; i < geneTrees.taxonIdToLabel.length; i++) {
            if (i > 0) System.out.print(", ");
            System.out.print(geneTrees.taxonIdToLabel[i]);
        }
        System.out.println();
        
        // Print STBipartitions with counts
        System.out.println("STBipartitions:");
        for (var entry : geneTrees.stBipartitions.entrySet()) {
            System.out.println("  " + entry.getKey().print(geneTrees.taxonIdToLabel) + " : " + entry.getValue());
        }
    }

    /**
     * Writes analysis results to the specified output file.
     * 
     * @param outputFilePath Path to the output file
     * @param content Content to write to the file
     * @throws IOException if there's an error writing to the file
     */
    private static void writeResults(String outputFilePath, String content) throws IOException {
        FileWriter writer = new FileWriter(outputFilePath);
        writer.write(content);
        writer.close();
    }
}
