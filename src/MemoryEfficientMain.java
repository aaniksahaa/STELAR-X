import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import preprocessing.MemoryEfficientGeneTrees;
import utils.Config;
import core.InferenceDP;  // Use original inference
import tree.STBipartition;  // Use original BitSet-based bipartitions
import tree.Tree;

/**
 * Memory-efficient main entry point that uses range-based processing during
 * gene tree preprocessing to reduce memory usage, but keeps the rest of the
 * pipeline (weight calculation, DP inference) exactly the same as the original.
 */
public class MemoryEfficientMain {

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

        // Parse command line arguments (same as original Main)
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
            } else if (args[i].equals("-e") && i + 1 < args.length) {
                expansionMethod = args[i + 1];
                i++;
            } else if (args[i].equals("-d") && i + 1 < args.length) {
                distanceMethod = args[i + 1];
                i++;
            } else if (args[i].equals("-v")) {
                verboseExpansion = true;
            } else if (args[i].equals("--no-expansion")) {
                disableExpansion = true;
            } else if (args[i].equals("-bs") && i + 1 < args.length) {
                branchSupport = args[i + 1];
                i++;
            } else if (args[i].equals("-lambda") && i + 1 < args.length) {
                lambda = Double.parseDouble(args[i + 1]);
                i++;
            }
        }

        // Validate required arguments
        if (inputFilePath == null) {
            System.err.println("Error: Input file (-i) is required");
            System.err.println("Usage: java MemoryEfficientMain -i input.tre [-o output.tre] [-m computation_mode] [other options]");
            return;
        }

        // Set defaults
        if (outputFilePath == null) {
            outputFilePath = "memory_efficient_output.tre";
        }
        if (computationMode == null) {
            computationMode = "CPU_PARALLEL";
        }
        if (expansionMethod == null) {
            expansionMethod = "NONE";  // Default to no expansion for simplicity
        }

        // Configure computation mode
        try {
            Config.COMPUTATION_MODE = Config.ComputationMode.valueOf(computationMode.toUpperCase());
        } catch (IllegalArgumentException e) {
            System.err.println("Invalid computation mode: " + computationMode);
            System.err.println("Valid modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
            return;
        }

        // Configure expansion (simplified - use defaults for most settings)
        if (disableExpansion || expansionMethod.equals("NONE")) {
            // Disable expansion by not setting any specific configuration
            System.out.println("Expansion disabled");
        }

        System.out.println("\n=== Memory-Efficient STELAR-MP ===");
        System.out.println("Input file: " + inputFilePath);
        System.out.println("Output file: " + outputFilePath);
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Expansion method: " + expansionMethod);
        System.out.println("Memory optimization: ENABLED (gene tree processing only)");
        System.out.println("===================================\n");

        try {
            long startTime = System.currentTimeMillis();

            // Phase 1: Memory-efficient gene trees preprocessing
            System.out.println("Phase 1: Loading gene trees with memory-efficient preprocessing...");
            long phase1Start = System.currentTimeMillis();
            
            MemoryEfficientGeneTrees memEffGeneTrees = new MemoryEfficientGeneTrees(inputFilePath);
            memEffGeneTrees.readTaxaNames();
            memEffGeneTrees.readGeneTrees(null);
            
            long phase1Time = System.currentTimeMillis() - phase1Start;
            System.out.println("Phase 1 completed in " + phase1Time + " ms");

            // Phase 2: Generate candidate bipartitions (same as original)
            System.out.println("\nPhase 2: Generating candidate bipartitions...");
            long phase2Start = System.currentTimeMillis();
            
            boolean useExpansion = !disableExpansion && !expansionMethod.equals("NONE");
            List<STBipartition> candidates = memEffGeneTrees.generateCandidateBipartitions(useExpansion);
            
            long phase2Time = System.currentTimeMillis() - phase2Start;
            System.out.println("Phase 2 completed in " + phase2Time + " ms");
            System.out.println("Generated " + candidates.size() + " candidate bipartitions");

            // Phase 3: Inference using original DP (no changes)
            System.out.println("\nPhase 3: Running inference with dynamic programming...");
            long phase3Start = System.currentTimeMillis();
            
            InferenceDP inference = new InferenceDP(memEffGeneTrees, candidates);
            double optimalScore = inference.solve();
            
            long phase3Time = System.currentTimeMillis() - phase3Start;
            System.out.println("Phase 3 completed in " + phase3Time + " ms");
            System.out.println("Optimal score: " + optimalScore);

            // Phase 4: Tree reconstruction (same as original)
            System.out.println("\nPhase 4: Reconstructing phylogenetic tree...");
            long phase4Start = System.currentTimeMillis();
            
            Tree resultTree = inference.reconstructTree();
            
            long phase4Time = System.currentTimeMillis() - phase4Start;
            System.out.println("Phase 4 completed in " + phase4Time + " ms");

            // Phase 5: Output
            System.out.println("\nPhase 5: Writing output...");
            String newickOutput = resultTree.getNewickFormat();
            
            try (FileWriter writer = new FileWriter(outputFilePath)) {
                writer.write(newickOutput);
            }

            long totalTime = System.currentTimeMillis() - startTime;
            
            // Summary
            System.out.println("\n=== Execution Summary ===");
            System.out.println("Phase 1 (Memory-efficient preprocessing): " + phase1Time + " ms");
            System.out.println("Phase 2 (Candidate generation): " + phase2Time + " ms");
            System.out.println("Phase 3 (Inference): " + phase3Time + " ms");
            System.out.println("Phase 4 (Tree reconstruction): " + phase4Time + " ms");
            System.out.println("Total execution time: " + totalTime + " ms");
            System.out.println("Optimal score: " + optimalScore);
            System.out.println("Output written to: " + outputFilePath);
            
            // Memory optimization summary
            System.out.println("\n=== Memory Optimization Summary ===");
            System.out.println("• Used range-based representation during gene tree processing");
            System.out.println("• Identified unique bipartitions before converting to BitSets");
            System.out.println("• Reduced memory usage during preprocessing phase");
            System.out.println("• Kept original inference pipeline unchanged for reliability");
            System.out.println("• Final unique bipartitions: " + memEffGeneTrees.stBipartitions.size());
            System.out.println("===================================");

        } catch (FileNotFoundException e) {
            System.err.println("Error: Input file not found: " + inputFilePath);
            System.err.println("Please check the file path and try again.");
        } catch (Exception e) {
            System.err.println("Error during execution: " + e.getMessage());
            e.printStackTrace();
        }
    }
}
