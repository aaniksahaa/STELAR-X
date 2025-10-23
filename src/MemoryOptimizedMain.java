import java.io.FileNotFoundException;
import java.util.List;

import core.MemoryEfficientInferenceDP;
import preprocessing.GeneTrees;
import tree.STBipartition;
import tree.Tree;
import utils.Config;

/**
 * Main class for memory-optimized STELAR implementation.
 * 
 * This version uses:
 * 1. Compact integer tuple representation for bipartitions
 * 2. Inverse index mapping for weight calculation
 * 3. Hash-based cluster representation for DP
 * 4. Memory-efficient GPU kernels (when available)
 * 
 * Memory usage is reduced from O(n²k) to O(nk) while maintaining
 * the same algorithmic complexity and accuracy.
 */
public class MemoryOptimizedMain {
    
    public static void main(String[] args) {
        if (args.length < 2) {
            System.err.println("Usage: java MemoryOptimizedMain <gene_trees_file> <output_file> [computation_mode]");
            System.err.println("Computation modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
            System.exit(1);
        }
        
        String geneTreesFile = args[0];
        String outputFile = args[1];
        
        // Set computation mode
        if (args.length >= 3) {
            try {
                Config.COMPUTATION_MODE = Config.ComputationMode.valueOf(args[2].toUpperCase());
            } catch (IllegalArgumentException e) {
                System.err.println("Invalid computation mode: " + args[2]);
                System.err.println("Valid modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
                System.exit(1);
            }
        } else {
            Config.COMPUTATION_MODE = Config.ComputationMode.CPU_PARALLEL; // Default
        }
        
        System.out.println("==== MEMORY-OPTIMIZED STELAR STARTED ====");
        System.out.println("Gene trees file: " + geneTreesFile);
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        
        try {
            long startTime = System.currentTimeMillis();
            
            System.out.println("==== MEMORY-OPTIMIZED STELAR EXECUTION STARTED ====");
            System.out.println("Input file: " + geneTreesFile);
            System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
            System.out.println("Java version: " + System.getProperty("java.version"));
            System.out.println("Available memory: " + (Runtime.getRuntime().maxMemory() / (1024*1024)) + " MB");
            
            // Step 1: Load and preprocess gene trees using memory-efficient approach
            System.out.println("\n=== STEP 1: GENE TREE PREPROCESSING ===");
            System.out.println("Loading gene trees from: " + geneTreesFile);
            
            GeneTrees geneTrees = new GeneTrees(geneTreesFile);
            System.out.println("GeneTrees object created successfully");
            
            System.out.println("Reading taxa names...");
            geneTrees.readTaxaNames();
            System.out.println("Taxa names read successfully");
            
            System.out.println("Reading gene trees...");
            geneTrees.readGeneTrees(null);
            System.out.println("Gene trees read successfully");
            
            long preprocessingTime = System.currentTimeMillis();
            System.out.println("Preprocessing completed in " + (preprocessingTime - startTime) + " ms");
            System.out.println("Taxa count: " + geneTrees.realTaxaCount);
            System.out.println("Gene trees: " + geneTrees.geneTrees.size());
            System.out.println("Unique bipartitions: " + geneTrees.stBipartitions.size());
            
            // Step 2: Generate candidate bipartitions
            System.out.println("\n=== STEP 2: CANDIDATE GENERATION ===");
            List<STBipartition> candidates = geneTrees.generateCandidateBipartitions();
            
            long candidateTime = System.currentTimeMillis();
            System.out.println("Candidate generation completed in " + (candidateTime - preprocessingTime) + " ms");
            System.out.println("Candidate bipartitions: " + candidates.size());
            
            // Debug: Check if we have candidates that span all taxa
            System.out.println("\nDEBUG: Analyzing candidate coverage...");
            int fullSpanCandidates = 0;
            for (STBipartition candidate : candidates) {
                int totalTaxa = candidate.cluster1.cardinality() + candidate.cluster2.cardinality();
                if (totalTaxa == geneTrees.realTaxaCount) {
                    fullSpanCandidates++;
                }
            }
            System.out.println("Candidates spanning all " + geneTrees.realTaxaCount + " taxa: " + fullSpanCandidates);
            
            if (fullSpanCandidates == 0) {
                System.err.println("ERROR: No candidates span all taxa! This will cause DP to fail.");
                System.err.println("Expected: At least 1 candidate with " + geneTrees.realTaxaCount + " total taxa");
                System.err.println("This suggests an issue in candidate generation.");
            }
            
            // Step 3: Memory-efficient DP inference
            System.out.println("\n=== STEP 3: MEMORY-EFFICIENT DP INFERENCE ===");
            System.out.println("Initializing memory-efficient DP with " + candidates.size() + " candidates...");
            
            System.out.println("About to create MemoryEfficientInferenceDP...");
            System.out.flush(); // Force output
            
            MemoryEfficientInferenceDP inference = null;
            try {
                inference = new MemoryEfficientInferenceDP(geneTrees, candidates);
                System.out.println("MemoryEfficientInferenceDP created successfully!");
            } catch (Exception e) {
                System.err.println("ERROR during MemoryEfficientInferenceDP creation:");
                System.err.println("  Error: " + e.getMessage());
                System.err.println("  Type: " + e.getClass().getSimpleName());
                e.printStackTrace();
                throw e;
            }
            
            System.out.println("Starting DP solve...");
            double optimalScore = inference.solve();
            System.out.println("DP solve completed with optimal score: " + optimalScore);
            
            long dpTime = System.currentTimeMillis();
            System.out.println("DP inference completed in " + (dpTime - candidateTime) + " ms");
            System.out.println("Optimal score: " + optimalScore);
            
            // Step 4: Tree reconstruction
            System.out.println("\n=== STEP 4: TREE RECONSTRUCTION ===");
            Tree resultTree = inference.reconstructTree();
            
            long reconstructionTime = System.currentTimeMillis();
            System.out.println("Tree reconstruction completed in " + (reconstructionTime - dpTime) + " ms");
            
            // Step 5: Output results
            System.out.println("\n=== STEP 5: OUTPUT ===");
            String newickTree = resultTree.getNewickFormat();
            System.out.println("Writing result tree to file...");
            
            // Write only the tree to the output file
            try (java.io.PrintWriter writer = new java.io.PrintWriter(outputFile)) {
                writer.println(newickTree);
            }
            
            System.out.println("Tree written to: " + outputFile);
            
            // Print final statistics
            long totalTime = System.currentTimeMillis() - startTime;
            System.out.println("\n==== MEMORY-OPTIMIZED STELAR COMPLETED ====");
            System.out.println("Total execution time: " + totalTime + " ms");
            System.out.println("Breakdown:");
            System.out.println("  Preprocessing: " + (preprocessingTime - startTime) + " ms");
            System.out.println("  Candidate generation: " + (candidateTime - preprocessingTime) + " ms");
            System.out.println("  DP inference: " + (dpTime - candidateTime) + " ms");
            System.out.println("  Tree reconstruction: " + (reconstructionTime - dpTime) + " ms");
            
            // Memory usage information
            Runtime runtime = Runtime.getRuntime();
            long usedMemory = runtime.totalMemory() - runtime.freeMemory();
            System.out.println("Memory usage: " + (usedMemory / (1024 * 1024)) + " MB");
            
            System.out.println("\nMemory optimization benefits:");
            System.out.println("- Reduced memory usage from O(n²k) to O(nk)");
            System.out.println("- Avoided bitset expansions during weight calculation and DP");
            System.out.println("- Used compact integer tuple representation throughout");
            System.out.println("- Employed hash-based cluster representation for DP memoization");
            
        } catch (FileNotFoundException e) {
            System.err.println("Error: Gene trees file not found: " + geneTreesFile);
            System.exit(1);
        } catch (Exception e) {
            System.err.println("Error during execution: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }
}
