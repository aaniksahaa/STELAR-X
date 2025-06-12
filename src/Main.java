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

        // Parse command line arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i") && i + 1 < args.length) {
                inputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            } else if (args[i].equals("-o") && i + 1 < args.length) {
                outputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            }
        }

        // Validate required arguments
        if (inputFilePath == null || outputFilePath == null) {
            System.out.println("Usage: java Main -i <input_file> -o <output_file>");
            System.exit(-1);
        }

        // Validate input file exists
        File inputFile = new File(inputFilePath);
        if (!inputFile.exists()) {
            System.err.println("Error: Input file '" + inputFilePath + "' does not exist.");
            System.exit(-1);
        }

        System.out.println("Input file: " + inputFilePath);
        System.out.println("Output file: " + outputFilePath);

        try {
            // Process gene trees and generate analysis results
            String analysisResult = processGeneTrees(inputFilePath);
            
            // Write results to output file
            writeResults(outputFilePath, analysisResult);
            
            System.out.println("Analysis complete! Results written to: " + outputFilePath);

        } catch (Exception e) {
            System.err.println("Error during analysis: " + e.getMessage());
            e.printStackTrace();
            System.exit(-1);
        }
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
        
        // Build analysis result with proper spacing from Config
        StringBuilder result = new StringBuilder();
        
        // Add n spaces before output
        for (int i = 0; i < Config.n; i++) {
            result.append(" ");
        }
        
        result.append("=== Gene Tree Analysis Results ===\n");
        
        for (int i = 0; i < Config.n; i++) {
            result.append(" ");
        }
        result.append("Number of gene trees: ").append(geneTreeCount).append("\n");
        
        for (int i = 0; i < Config.n; i++) {
            result.append(" ");
        }
        result.append("Number of taxa: ").append(taxaCount).append("\n");
        
        for (int i = 0; i < Config.n; i++) {
            result.append(" ");
        }
        result.append("Unique tripartitions: ").append(uniquePartitions).append("\n");
        
        for (int i = 0; i < Config.n; i++) {
            result.append(" ");
        }
        result.append("Unique STBipartitions: ").append(uniqueSTBipartitions).append("\n");
        
        // Add taxa names
        for (int i = 0; i < Config.n; i++) {
            result.append(" ");
        }
        result.append("Taxa names: ");
        for (int i = 0; i < geneTrees.taxonIdToLabel.length; i++) {
            if (i > 0) result.append(", ");
            result.append(geneTrees.taxonIdToLabel[i]);
        }
        result.append("\n");
        
        // Add STBipartitions with counts
        for (int i = 0; i < Config.n; i++) {
            result.append(" ");
        }
        result.append("STBipartitions:\n");
        for (var entry : geneTrees.stBipartitions.entrySet()) {
            for (int i = 0; i < Config.n; i++) {
                result.append(" ");
            }
            result.append("  ").append(entry.getKey().print(geneTrees.taxonIdToLabel))
                  .append(" : ").append(entry.getValue()).append("\n");
        }
        
        // Test InferenceDP algorithm
        System.out.println("Testing InferenceDP algorithm...");
        List<STBipartition> candidates = new ArrayList<>(geneTrees.stBipartitions.keySet());
        
        if (!candidates.isEmpty()) {
            InferenceDP dp = new InferenceDP(geneTrees, candidates);
            double maxScore = dp.solve();
            
            for (int i = 0; i < Config.n; i++) {
                result.append(" ");
            }
            result.append("DP Algorithm Results:\n");
            for (int i = 0; i < Config.n; i++) {
                result.append(" ");
            }
            result.append("  Maximum score: ").append(maxScore).append("\n");
            
            Tree reconstructedTree = dp.reconstructTree();
            if (reconstructedTree != null && reconstructedTree.root != null) {
                for (int i = 0; i < Config.n; i++) {
                    result.append(" ");
                }
                result.append("  Reconstructed tree: ").append(reconstructedTree.getNewickFormat()).append("\n");
            }
        }
        
        return result.toString();
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
