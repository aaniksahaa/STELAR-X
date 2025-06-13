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
import utils.BitSet;

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
            }
        }

        // Validate required arguments
        if (inputFilePath == null || outputFilePath == null) {
            System.out.println("Usage: java Main -i <input_file> -o <output_file> [-m <computation_mode>]");
            System.out.println("Computation modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
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

        System.out.println("Input file: " + inputFilePath);
        System.out.println("Output file: " + outputFilePath);
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);

        long startTime = System.nanoTime();

        // Process gene trees
        GeneTrees geneTrees = new GeneTrees(inputFilePath);
        geneTrees.readTaxaNames(); // Ensure taxaMap is initialized
        geneTrees.readGeneTrees(null);

        // Generate candidate bipartitions
        List<STBipartition> candidates = geneTrees.generateCandidateBipartitions();

        // Run inference
        InferenceDP inference = new InferenceDP(geneTrees, candidates);
        double score = inference.solve();
        Tree resultTree = inference.reconstructTree();

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
        int uniquePartitions = geneTrees.triPartitions.size();
        int uniqueSTBipartitions = geneTrees.stBipartitions.size();
        
        System.out.println("Processing complete:");
        System.out.println("  - Gene trees processed: " + geneTreeCount);
        System.out.println("  - Taxa found: " + taxaCount);
        System.out.println("  - Unique tripartitions: " + uniquePartitions);
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
