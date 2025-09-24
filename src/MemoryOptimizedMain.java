import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import preprocessing.MemoryOptimizedGeneTrees;
import core.MemoryOptimizedInferenceDP;
import tree.RangeSTBipartition;
import tree.Tree;
import utils.Config;

/**
 * Memory-optimized main entry point for phylogeny reconstruction.
 * 
 * This implementation uses range-based bipartition representation to dramatically
 * reduce memory usage from O(n²k) to O(nk) where n=taxa, k=gene trees.
 * 
 * Key optimizations:
 * 1. Range-based bipartition representation instead of BitSets
 * 2. Hash-based equality checking for cross-tree bipartition comparison
 * 3. Efficient GPU computation using array traversal instead of bitset operations
 * 4. Proper multithreading with reduced memory overhead
 */
public class MemoryOptimizedMain {

    public static void main(String[] args) throws IOException {
        
        String inputFilePath = null;
        String outputFilePath = null;
        String computationMode = null;
        double lambda = 0.5;
        boolean enableProfiling = false;
        
        // Parse command line arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i") && i + 1 < args.length) {
                inputFilePath = args[i + 1];
                i++;
            } else if (args[i].equals("-o") && i + 1 < args.length) {
                outputFilePath = args[i + 1];
                i++;
            } else if (args[i].equals("-m") && i + 1 < args.length) {
                computationMode = args[i + 1];
                i++;
            } else if (args[i].equals("-lambda") && i + 1 < args.length) {
                lambda = Double.parseDouble(args[i + 1]);
                i++;
            } else if (args[i].equals("-profile")) {
                enableProfiling = true;
            } else if (args[i].equals("-h") || args[i].equals("--help")) {
                printUsage();
                return;
            }
        }

        // Validate required arguments
        if (inputFilePath == null) {
            System.err.println("Error: Input file path is required (-i)");
            printUsage();
            return;
        }

        if (outputFilePath == null) {
            outputFilePath = "memory_optimized_output.tre";
            System.out.println("Using default output file: " + outputFilePath);
        }

        // Set computation mode
        if (computationMode == null) {
            computationMode = "CPU_PARALLEL";
        }
        
        try {
            Config.COMPUTATION_MODE = Config.ComputationMode.valueOf(computationMode.toUpperCase());
        } catch (IllegalArgumentException e) {
            System.err.println("Invalid computation mode: " + computationMode);
            System.err.println("Valid modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
            return;
        }

        System.out.println("\n=== Memory-Optimized Phylogeny Reconstruction ===");
        System.out.println("Input file: " + inputFilePath);
        System.out.println("Output file: " + outputFilePath);
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Lambda parameter: " + lambda);
        System.out.println("Memory optimization: ENABLED");
        System.out.println("===============================================\n");

        try {
            long startTime = System.currentTimeMillis();
            
            // Phase 1: Memory-optimized gene tree preprocessing
            System.out.println("Phase 1: Loading and preprocessing gene trees with memory optimization...");
            long phase1Start = System.currentTimeMillis();
            
            MemoryOptimizedGeneTrees memOptGeneTrees = new MemoryOptimizedGeneTrees(inputFilePath);
            memOptGeneTrees.readTaxaNames();
            memOptGeneTrees.readGeneTrees(null);
            
            long phase1Time = System.currentTimeMillis() - phase1Start;
            System.out.println("Phase 1 completed in " + phase1Time + " ms");
            
            if (enableProfiling) {
                printMemoryUsage("After Phase 1");
            }

            // Phase 2: Generate range-based candidate bipartitions
            System.out.println("\nPhase 2: Generating range-based candidate bipartitions...");
            long phase2Start = System.currentTimeMillis();
            
            List<RangeSTBipartition> rangeCandidates = memOptGeneTrees.generateRangeCandidateBipartitions();
            
            long phase2Time = System.currentTimeMillis() - phase2Start;
            System.out.println("Phase 2 completed in " + phase2Time + " ms");
            System.out.println("Generated " + rangeCandidates.size() + " range-based candidates");
            
            if (enableProfiling) {
                printMemoryUsage("After Phase 2");
            }

            // Phase 3: Memory-optimized inference with dynamic programming
            System.out.println("\nPhase 3: Memory-optimized inference with dynamic programming...");
            long phase3Start = System.currentTimeMillis();
            
            MemoryOptimizedInferenceDP memOptInference = 
                new MemoryOptimizedInferenceDP(memOptGeneTrees, rangeCandidates);
            
            double optimalScore = memOptInference.solve();
            
            long phase3Time = System.currentTimeMillis() - phase3Start;
            System.out.println("Phase 3 completed in " + phase3Time + " ms");
            System.out.println("Optimal score: " + optimalScore);
            
            if (enableProfiling) {
                printMemoryUsage("After Phase 3");
                memOptInference.printMemoryStats();
            }

            // Phase 4: Tree reconstruction
            System.out.println("\nPhase 4: Reconstructing phylogenetic tree...");
            long phase4Start = System.currentTimeMillis();
            
            Tree resultTree = memOptInference.reconstructTree();
            
            long phase4Time = System.currentTimeMillis() - phase4Start;
            System.out.println("Phase 4 completed in " + phase4Time + " ms");

            // Phase 5: Output tree in Newick format
            System.out.println("\nPhase 5: Writing output tree...");
            String newickOutput = resultTree.getNewickFormat();
            
            try (FileWriter writer = new FileWriter(outputFilePath)) {
                writer.write(newickOutput);
            }
            
            long totalTime = System.currentTimeMillis() - startTime;
            
            // Final statistics
            System.out.println("\n=== Execution Summary ===");
            System.out.println("Phase 1 (Preprocessing): " + phase1Time + " ms");
            System.out.println("Phase 2 (Candidate generation): " + phase2Time + " ms");
            System.out.println("Phase 3 (Inference): " + phase3Time + " ms");
            System.out.println("Phase 4 (Tree reconstruction): " + phase4Time + " ms");
            System.out.println("Total execution time: " + totalTime + " ms");
            System.out.println("Optimal score: " + optimalScore);
            System.out.println("Output written to: " + outputFilePath);
            System.out.println("=========================");

            // Final memory usage
            if (enableProfiling) {
                printMemoryUsage("Final");
                System.out.println("\n=== Memory Optimization Benefits ===");
                System.out.println("Traditional representation: O(n²k) memory");
                System.out.println("Range-based representation: O(nk) memory");
                System.out.println("Estimated memory reduction: " + 
                    (memOptGeneTrees.realTaxaCount > 0 ? 
                     (double)memOptGeneTrees.realTaxaCount + "x" : "N/A"));
                System.out.println("====================================");
            }

        } catch (FileNotFoundException e) {
            System.err.println("Error: Input file not found: " + inputFilePath);
            System.err.println("Please check the file path and try again.");
        } catch (Exception e) {
            System.err.println("Error during execution: " + e.getMessage());
            e.printStackTrace();
            
            if (enableProfiling) {
                printMemoryUsage("Error state");
            }
        }
    }
    
    /**
     * Print current memory usage statistics.
     */
    private static void printMemoryUsage(String phase) {
        Runtime runtime = Runtime.getRuntime();
        long totalMemory = runtime.totalMemory();
        long freeMemory = runtime.freeMemory();
        long usedMemory = totalMemory - freeMemory;
        long maxMemory = runtime.maxMemory();
        
        System.out.println("\n--- Memory Usage (" + phase + ") ---");
        System.out.println("Used memory: " + (usedMemory / 1024 / 1024) + " MB");
        System.out.println("Total memory: " + (totalMemory / 1024 / 1024) + " MB");
        System.out.println("Max memory: " + (maxMemory / 1024 / 1024) + " MB");
        System.out.println("Free memory: " + (freeMemory / 1024 / 1024) + " MB");
        System.out.println("Memory utilization: " + 
            String.format("%.1f%%", (double)usedMemory / maxMemory * 100));
        System.out.println("---------------------------\n");
    }

    /**
     * Print usage information.
     */
    private static void printUsage() {
        System.out.println("Memory-Optimized Phylogeny Reconstruction");
        System.out.println("Usage: java MemoryOptimizedMain [options]");
        System.out.println();
        System.out.println("Required arguments:");
        System.out.println("  -i <file>        Input file path (gene trees in Newick format)");
        System.out.println();
        System.out.println("Optional arguments:");
        System.out.println("  -o <file>        Output file path (default: memory_optimized_output.tre)");
        System.out.println("  -m <mode>        Computation mode: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
        System.out.println("                   (default: CPU_PARALLEL)");
        System.out.println("  -lambda <value>  Lambda parameter for scoring (default: 0.5)");
        System.out.println("  -profile         Enable detailed memory and performance profiling");
        System.out.println("  -h, --help       Show this help message");
        System.out.println();
        System.out.println("Memory Optimization Features:");
        System.out.println("  • Range-based bipartition representation (O(1) space per bipartition)");
        System.out.println("  • Hash-based equality checking for efficient uniqueness detection");
        System.out.println("  • GPU-optimized array traversal instead of bitset operations");
        System.out.println("  • Reduced memory footprint from O(n²k) to O(nk)");
        System.out.println();
        System.out.println("Example:");
        System.out.println("  java MemoryOptimizedMain -i input_trees.tre -o result.tre -m CPU_PARALLEL -profile");
    }
}
