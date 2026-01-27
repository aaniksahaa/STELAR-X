package core;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;

import com.sun.jna.Memory;

import preprocessing.GeneTrees;
import tree.*;
import utils.Config;
import utils.Threading;
import core.WeightCalculator.WeightCalcLib;

/**
 * Memory-optimized weight calculator that avoids BitSet expansion.
 * 
 * This implementation uses:
 * 1. Range-based bipartition representations from MemoryEfficientBipartitionManager
 * 2. Inverse index arrays for O(min(|A|, |B|)) intersection calculations
 * 3. Direct range processing without BitSet conversion
 * 
 * Memory usage: O(nk) instead of O(n²k) compared to traditional BitSet approach.
 * 
 * Supports all computation modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL
 */
public class MemoryOptimizedWeightCalculator {
    
    private final GeneTrees geneTrees;
    private final InverseIndexManager inverseIndexManager;
    private final MemoryEfficientBipartitionManager bipartitionManager;
    private final Map<Object, List<RangeBipartition>> hashToBipartitions;
    
    // Statistics for performance monitoring
    private long totalScoreCalculations = 0;
    private long totalIntersectionCalculations = 0;
    private long totalProcessingTime = 0;
    
    public MemoryOptimizedWeightCalculator(GeneTrees geneTrees) {
        System.out.println("==== INITIALIZING MEMORY-OPTIMIZED WEIGHT CALCULATOR ====");
        
        this.geneTrees = geneTrees;
        
        System.out.println("Creating MemoryEfficientBipartitionManager...");
        this.bipartitionManager = new MemoryEfficientBipartitionManager(
            geneTrees.geneTrees, geneTrees.realTaxaCount);
        
        System.out.println("Processing gene trees to extract range bipartitions...");
        bipartitionManager.processGeneTreesParallel();
        
        System.out.println("Creating InverseIndexManager...");
        this.inverseIndexManager = new InverseIndexManager(
            geneTrees.geneTrees, geneTrees.realTaxaCount);
        
        this.hashToBipartitions = bipartitionManager.getHashToBipartitions();
        
        System.out.println("Memory-optimized weight calculator initialized");
        System.out.println("Range bipartition groups: " + hashToBipartitions.size());
        System.out.println("Total range bipartitions: " + 
                         hashToBipartitions.values().stream().mapToInt(List::size).sum());
        System.out.println("==== MEMORY-OPTIMIZED WEIGHT CALCULATOR READY ====");
    }
    
    /**
     * Calculate weights for candidate bipartitions using memory-optimized approach.
     */
    public Map<RangeBipartition, Double> calculateWeights(List<RangeBipartition> candidates) {
        System.out.println("==== MEMORY-OPTIMIZED WEIGHT CALCULATION STARTED ====");
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Number of candidates: " + candidates.size());
        System.out.println("Gene tree bipartition groups: " + hashToBipartitions.size());
        
        long startTime = System.currentTimeMillis();
        
        Map<RangeBipartition, Double> result;
        switch (Config.COMPUTATION_MODE) {
            case CPU_SINGLE:
                System.out.println("Using memory-optimized CPU_SINGLE mode");
                result = calculateWeightsSingleThread(candidates);
                break;
            case CPU_PARALLEL:
                System.out.println("Using memory-optimized CPU_PARALLEL mode");
                result = calculateWeightsMultiThread(candidates);
                break;
            case GPU_PARALLEL:
                System.out.println("Using memory-optimized GPU_PARALLEL mode");
                result = calculateWeightsCompactGPU(candidates);
                break;
            default:
                throw new IllegalStateException("Unknown computation mode: " + Config.COMPUTATION_MODE);
        }
        
        long endTime = System.currentTimeMillis();
        totalProcessingTime = endTime - startTime;
        
        System.out.println("==== MEMORY-OPTIMIZED WEIGHT CALCULATION COMPLETED ====");
        System.out.println("Processing time: " + totalProcessingTime + " ms");
        System.out.println("Total score calculations: " + totalScoreCalculations);
        System.out.println("Total intersection calculations: " + totalIntersectionCalculations);
        System.out.println("Average intersections per score: " + 
                         (totalScoreCalculations > 0 ? totalIntersectionCalculations / (double) totalScoreCalculations : 0));
        
        // Print inverse index statistics
        System.out.println("\n" + inverseIndexManager.getStatistics());
        
        return result;
    }
    
    /**
     * Single-threaded weight calculation using range-based processing.
     */
    private Map<RangeBipartition, Double> calculateWeightsSingleThread(List<RangeBipartition> candidates) {
        System.out.println("Starting single-threaded range-based weight calculation...");
        
        Map<RangeBipartition, Double> weights = new HashMap<>();
        
        // Use direct calculation for candidates (small number) and range-based for gene trees (large number)
        System.out.println("Using hybrid approach: direct calculation for candidates, range-based for gene trees");
        
        int processedCandidates = 0;
        
        for (RangeBipartition candidate : candidates) {
            double totalScore = 0.0;
            
            // Use traditional BitSet calculation for candidates (simple and direct)
            totalScore = calculateRangeBasedScore(candidate);
            
            weights.put(candidate, totalScore);
            processedCandidates++;
            
            // Log progress for large datasets
            if (processedCandidates % 1000 == 0 || processedCandidates == candidates.size()) {
                System.out.println("Processed " + processedCandidates + "/" + candidates.size() + " candidates");
            }
        }
        
        System.out.println("Single-threaded calculation completed");
        return weights;
    }
    
    /**
     * Multi-threaded weight calculation using range-based processing.
     */
    private Map<RangeBipartition, Double> calculateWeightsMultiThread(List<RangeBipartition> candidates) {
        System.out.println("Starting multi-threaded range-based weight calculation...");
        
        Map<RangeBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        Threading.startThreading(numThreads);
        
        // Use hybrid approach: direct calculation for candidates, range-based for gene trees
        System.out.println("Using hybrid approach for multi-threaded calculation");
        
        // Calculate optimal number of threads to avoid invalid ranges
        int chunkSize = Math.max(1, (candidates.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (candidates.size() + chunkSize - 1) / chunkSize);
        CountDownLatch latch = new CountDownLatch(actualThreads);
        
        System.out.println("Using " + actualThreads + " threads with chunk size " + chunkSize);
        
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, candidates.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    // Validate range before processing
                    if (startIdx >= endIdx || startIdx >= candidates.size()) {
                        System.out.println("Thread " + threadId + " skipped - invalid range [" + startIdx + ", " + endIdx + ")");
                        return;
                    }
                    
                    System.out.println("Thread " + threadId + " processing candidates " + startIdx + " to " + (endIdx - 1));
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        RangeBipartition candidate = candidates.get(j);
                        double totalScore = 0.0;
                        
                        // Use traditional BitSet calculation (simple and direct)
                        totalScore = calculateRangeBasedScore(candidate);
                        
                        weights.put(candidate, totalScore);
                    }
                    
                    System.out.println("Thread " + threadId + " completed processing " + (endIdx - startIdx) + " candidates");
                    
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Weight calculation was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        System.out.println("Multi-threaded calculation completed");
        return weights;
    }
    
    
    /**
     * Calculate score between two range bipartitions using inverse index.
     * 
     * ENHANCED: Properly handles trees with different taxa using sentinel values.
     * - InverseIndexManager automatically handles -1 sentinel values
     * - Only counts intersections for taxa present in both trees
     * 
     * Implements the same scoring formula as original calculateScore method.
     */
    private double calculateRangeScore(RangeBipartition range1, RangeBipartition range2) {
        // Calculate four intersection sizes: AA, AB, BA, BB
        // InverseIndexManager handles sentinel values automatically
        int aa = inverseIndexManager.getRangeIntersectionSize(
            range1.geneTreeIndex, range1.leftStart, range1.leftEnd,
            range2.geneTreeIndex, range2.leftStart, range2.leftEnd);
            
        int bb = inverseIndexManager.getRangeIntersectionSize(
            range1.geneTreeIndex, range1.rightStart, range1.rightEnd,
            range2.geneTreeIndex, range2.rightStart, range2.rightEnd);
            
        int ab = inverseIndexManager.getRangeIntersectionSize(
            range1.geneTreeIndex, range1.leftStart, range1.leftEnd,
            range2.geneTreeIndex, range2.rightStart, range2.rightEnd);
            
        int ba = inverseIndexManager.getRangeIntersectionSize(
            range1.geneTreeIndex, range1.rightStart, range1.rightEnd,
            range2.geneTreeIndex, range2.leftStart, range2.leftEnd);
        
        totalIntersectionCalculations += 4;
        
        // Apply same scoring formula as original implementation
        double score1 = 0;
        if (aa + bb >= 2) {
            score1 = aa * bb * (aa + bb - 2) / 2.0;
        }
        
        double score2 = 0;
        if (ab + ba >= 2) {
            score2 = ab * ba * (ab + ba - 2) / 2.0;
        }
        
        return score1 + score2;
    }
    
    /**
     * Calculate score using range-based approach.
     * Uses the memory-optimized range intersection method.
     */
    private double calculateRangeBasedScore(RangeBipartition candidate) {
        double totalScore = 0.0;
        
        // Use the existing gene tree RangeBipartitions from GeneTrees
        for (Map.Entry<RangeBipartition, Integer> entry : geneTrees.rangeBipartitions.entrySet()) {
            RangeBipartition geneTreeRange = entry.getKey();
            int frequency = entry.getValue();
            
            double score = calculateRangeScore(candidate, geneTreeRange);
            totalScore += score * frequency;
            
            totalScoreCalculations++;
        }
        
        return totalScore;
    }
    
    
    // Removed expensive buildCandidateRangeMapping and findMatchingRange methods
    // Using direct calculation approach instead
    
    // Removed unused utility methods:
    // - bitSetToTaxonSet: No longer needed with range-based approach
    // - getRangeTaxonSet: Replaced by InverseIndexManager functionality
    
    /**
     * Memory-optimized GPU weight calculation using compact range representations.
     * This replaces the traditional BitSet-based GPU approach with O(nk) memory usage.
     */
    private Map<RangeBipartition, Double> calculateWeightsCompactGPU(List<RangeBipartition> candidates) {
        System.out.println("==== STARTING COMPACT GPU WEIGHT CALCULATION ====");
        
        try {
            // Convert candidates and gene tree bipartitions to compact representations
            System.out.println("Converting bipartitions to compact range representations...");
            
            // Use hybrid GPU approach: BitSet candidates vs Range gene trees
            System.out.println("Using hybrid GPU approach: BitSet candidates vs compact range gene trees");
            
            // Extract gene tree ranges and frequencies
            List<RangeBipartition> geneTreeRanges = new ArrayList<>();
            List<Integer> frequencies = new ArrayList<>();
            
            for (Map.Entry<Object, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
                List<RangeBipartition> ranges = entry.getValue();
                if (!ranges.isEmpty()) {
                    geneTreeRanges.add(ranges.get(0)); // Representative
                    frequencies.add(ranges.size());    // Frequency
                }
            }
            
            System.out.println("Gene tree ranges: " + geneTreeRanges.size());
            System.out.println("Candidates: " + candidates.size());
            
            // Check if CUDA library is available
            try {
                WeightCalcLib.INSTANCE.toString(); // Test if library loads
                System.out.println("CUDA library loaded successfully");
            } catch (UnsatisfiedLinkError e) {
                System.err.println("CUDA library not found: " + e.getMessage());
                System.out.println("Falling back to CPU calculation...");
                return calculateWeightsMultiThread(candidates);
            }
            
            // For candidates, we'll use a simpler approach:
            // Since candidates are small in number, use traditional BitSet GPU calculation
            // But use compact ranges for gene trees (which are large in number)
            System.out.println("Using hybrid approach: BitSet candidates vs compact gene tree ranges");
            
            return calculateWeightsHybridGPU(candidates, geneTreeRanges, frequencies);
            
        } catch (Exception e) {
            System.err.println("GPU calculation failed: " + e.getMessage());
            e.printStackTrace();
            System.out.println("Falling back to CPU calculation...");
            return calculateWeightsMultiThread(candidates);
        }
    }
    
    /**
     * Pure range-based GPU calculation using compact bipartitions.
     * This completely eliminates BitSet usage and uses the compact GPU kernel.
     */
    private Map<RangeBipartition, Double> calculateWeightsHybridGPU(
            List<RangeBipartition> candidates,
            List<RangeBipartition> geneTreeRanges,
            List<Integer> frequencies) {
        
        System.out.println("==== STARTING PURE RANGE-BASED GPU CALCULATION ====");
        System.out.println("Range candidates: " + candidates.size());
        System.out.println("Compact gene tree ranges: " + geneTreeRanges.size());
        
        try {
            int numCandidates = candidates.size();
            int numGeneTreeBips = geneTreeRanges.size();
            int numTrees = inverseIndexManager.getNumTrees();
            int numTaxa = inverseIndexManager.getNumTaxa();
            
            System.out.println("Pure Range GPU Parameters:");
            System.out.println("  Candidates: " + numCandidates);
            System.out.println("  Gene tree bipartitions: " + numGeneTreeBips);
            System.out.println("  Trees: " + numTrees);
            System.out.println("  Taxa: " + numTaxa);
            
            // Convert candidates to GPU compact format using contiguous memory allocation
            WeightCalcLib.CompactBipartition candidateTemplate = new WeightCalcLib.CompactBipartition();
            WeightCalcLib.CompactBipartition[] candidateArray = 
                (WeightCalcLib.CompactBipartition[]) candidateTemplate.toArray(numCandidates);
            
            for (int i = 0; i < numCandidates; i++) {
                RangeBipartition candidate = candidates.get(i);
                candidateArray[i].geneTreeIndex = candidate.geneTreeIndex;
                candidateArray[i].leftStart = candidate.leftStart;
                candidateArray[i].leftEnd = candidate.leftEnd;
                candidateArray[i].rightStart = candidate.rightStart;
                candidateArray[i].rightEnd = candidate.rightEnd;
            }
            
            // Convert gene tree ranges to GPU compact format using contiguous memory allocation
            WeightCalcLib.CompactBipartition geneTreeTemplate = new WeightCalcLib.CompactBipartition();
            WeightCalcLib.CompactBipartition[] geneTreeArray = 
                (WeightCalcLib.CompactBipartition[]) geneTreeTemplate.toArray(numGeneTreeBips);
            
            for (int i = 0; i < numGeneTreeBips; i++) {
                RangeBipartition range = geneTreeRanges.get(i);
                geneTreeArray[i].geneTreeIndex = range.geneTreeIndex;
                geneTreeArray[i].leftStart = range.leftStart;
                geneTreeArray[i].leftEnd = range.leftEnd;
                geneTreeArray[i].rightStart = range.rightStart;
                geneTreeArray[i].rightEnd = range.rightEnd;
            }
            
            // Prepare frequency array
            int[] frequencyArray = frequencies.stream().mapToInt(Integer::intValue).toArray();
            
            // Prepare result array
            double[] weights = new double[numCandidates];
            
            // Flatten inverse index and ordering arrays for GPU
            Memory inverseIndexMemory = flattenInverseIndex();
            Memory orderingMemory = flattenOrderings();
            
            // Launch compact GPU kernel
            System.out.println("Launching pure range-based GPU kernel with inverse index arrays...");
            System.out.println("IMPORTANT: GPU kernel must handle -1 sentinel values in inverse index");
            System.out.println("  - inverseIndex[tree][taxon] == -1 means taxon not present in tree");
            System.out.println("  - Only count intersections for taxa present in both trees");
            
            WeightCalcLib.INSTANCE.launchCompactWeightCalculation(
                candidateArray,
                geneTreeArray,
                frequencyArray,
                weights,
                inverseIndexMemory,
                orderingMemory,
                numCandidates,
                numGeneTreeBips,
                numTrees,
                numTaxa
            );
            
            System.out.println("==== PURE RANGE-BASED GPU KERNEL COMPLETED ====");
            
            // Convert results back to Java Map
            Map<RangeBipartition, Double> result = new HashMap<>();
            for (int i = 0; i < numCandidates; i++) {
                result.put(candidates.get(i), weights[i]);
            }
            
            return result;
            
        } catch (Exception e) {
            System.err.println("Pure range-based GPU calculation failed: " + e.getMessage());
            e.printStackTrace();
            System.out.println("Falling back to CPU calculation...");
            return calculateWeightsMultiThread(candidates);
        }
    }
    
    
    /**
     * Flatten inverse index array for GPU transfer.
     * 
     * CRITICAL: Handles sentinel values for trees with different taxa.
     * - GPU kernel MUST check for -1 values before using positions
     * - Format: memory[tree * numTaxa + taxon] = position (or -1 if taxon not in tree)
     * - Memory layout: [tree0_taxon0, tree0_taxon1, ..., tree1_taxon0, tree1_taxon1, ...]
     */
    @SuppressWarnings("resource") // Memory is managed by JNA and GPU kernel
    private Memory flattenInverseIndex() {
        int[][] inverseIndex = inverseIndexManager.getInverseIndex();
        int numTrees = inverseIndex.length;
        int numTaxa = inverseIndex[0].length;
        
        System.out.println("==== FLATTENING INVERSE INDEX WITH SENTINEL SUPPORT ====");
        System.out.println("Inverse index dimensions: " + numTrees + " trees x " + numTaxa + " taxa");
        System.out.println("Sentinel value: -1 (indicates taxon not present in tree)");
        
        Memory memory = new Memory((long) numTrees * numTaxa * 4); // 4 bytes per int
        
        int totalSentinels = 0;
        int totalValidPositions = 0;
        
        for (int tree = 0; tree < numTrees; tree++) {
            // Validate array dimensions
            if (inverseIndex[tree].length != numTaxa) {
                System.err.println("ERROR: Inverse index tree " + tree + " has " + 
                                 inverseIndex[tree].length + " taxa, expected " + numTaxa);
                throw new RuntimeException("Inconsistent inverse index dimensions");
            }
            
            int treeSentinels = 0;
            
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                long offset = ((long) tree * numTaxa + taxon) * 4;
                int position = inverseIndex[tree][taxon];
                
                memory.setInt(offset, position);
                
                // Count sentinels for statistics
                if (position == -1) {
                    treeSentinels++;
                    totalSentinels++;
                } else {
                    totalValidPositions++;
                }
            }
            
            // Log progress for large datasets (with sentinel statistics)
            if (tree % 100 == 0 || tree == numTrees - 1) {
                System.out.println("Flattened tree " + (tree + 1) + "/" + numTrees + 
                                 " (" + treeSentinels + " sentinels, " + 
                                 (numTaxa - treeSentinels) + " valid positions)");
            }
        }
        
        System.out.println("Inverse index flattening statistics:");
        System.out.println("  Total positions: " + (numTrees * numTaxa));
        System.out.println("  Valid positions: " + totalValidPositions);
        System.out.println("  Sentinel positions: " + totalSentinels);
        System.out.println("  Taxa coverage: " + String.format("%.2f%%", 
                         100.0 * totalValidPositions / (numTrees * numTaxa)));
        System.out.println("==== INVERSE INDEX FLATTENING COMPLETED ====");
        
        // Note: Memory object will be automatically freed by JNA when no longer referenced
        // GPU kernel is responsible for copying data before this method returns
        return memory;
    }
    
    /**
     * Flatten ordering arrays for GPU transfer.
     */
    @SuppressWarnings("resource") // Memory is managed by JNA and GPU kernel
    private Memory flattenOrderings() {
        int[][] orderings = inverseIndexManager.getGeneTreeOrderings();
        int numTrees = orderings.length;
        
        System.out.println("==== FLATTENING ORDERINGS WITH VARIABLE TREE SIZES ====");
        System.out.println("Number of trees: " + numTrees);
        
        // Find the maximum number of taxa across all trees
        int maxNumTaxa = 0;
        int minNumTaxa = Integer.MAX_VALUE;
        int totalTaxa = 0;
        
        for (int tree = 0; tree < numTrees; tree++) {
            if (orderings[tree] == null) {
                System.err.println("ERROR: Tree " + tree + " has null ordering array");
                throw new RuntimeException("Null ordering array for tree " + tree);
            }
            
            int actualLength = orderings[tree].length;
            maxNumTaxa = Math.max(maxNumTaxa, actualLength);
            minNumTaxa = Math.min(minNumTaxa, actualLength);
            totalTaxa += actualLength;
            
            // Log first few trees and some statistics
            if (tree < 5) {
                System.out.println("Tree " + tree + " length: " + actualLength);
            }
        }
        
        System.out.println("Tree size statistics:");
        System.out.println("  Min taxa per tree: " + minNumTaxa);
        System.out.println("  Max taxa per tree: " + maxNumTaxa);
        System.out.println("  Average taxa per tree: " + (totalTaxa / (double) numTrees));
        System.out.println("Using padded size: " + maxNumTaxa + " taxa per tree");
        
        // Use the maximum size and pad shorter trees with -1 (sentinel value)
        int paddedNumTaxa = maxNumTaxa;
        Memory memory = new Memory((long) numTrees * paddedNumTaxa * 4); // 4 bytes per int
        System.out.println("Allocated memory for " + numTrees + " x " + paddedNumTaxa + " = " + 
                         (numTrees * paddedNumTaxa) + " integers");
        
        int paddedPositions = 0;
        
        for (int tree = 0; tree < numTrees; tree++) {
            int actualTreeSize = orderings[tree].length;
            
            // Copy actual taxa positions
            for (int pos = 0; pos < actualTreeSize; pos++) {
                long offset = ((long) tree * paddedNumTaxa + pos) * 4;
                memory.setInt(offset, orderings[tree][pos]);
            }
            
            // Pad remaining positions with -1 (sentinel value indicating "no taxon")
            for (int pos = actualTreeSize; pos < paddedNumTaxa; pos++) {
                long offset = ((long) tree * paddedNumTaxa + pos) * 4;
                memory.setInt(offset, -1); // Sentinel value
                paddedPositions++;
            }
            
            // // Log progress for large datasets
            // if (tree % 100 == 0 || tree == numTrees - 1) {
            //     System.out.println("Processed " + (tree + 1) + "/" + numTrees + 
            //                      " trees (tree " + tree + " has " + actualTreeSize + " taxa)");
            // }
        }
        
        System.out.println("Padding statistics:");
        System.out.println("  Total padded positions: " + paddedPositions);
        System.out.println("  Padding overhead: " + (paddedPositions / (double)(numTrees * paddedNumTaxa) * 100) + "%");
        System.out.println("==== FLATTEN ORDERINGS COMPLETED ====");
        
        // Note: Memory object will be automatically freed by JNA when no longer referenced
        // GPU kernel is responsible for copying data before this method returns
        return memory;
    }
    
    // Removed expensive buildCandidateRangeMapping and findMatchingRange methods
    // Using direct calculation approach instead
}
