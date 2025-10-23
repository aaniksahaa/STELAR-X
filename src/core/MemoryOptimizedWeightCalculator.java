package core;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;

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
    public Map<STBipartition, Double> calculateWeights(List<STBipartition> candidates) {
        System.out.println("==== MEMORY-OPTIMIZED WEIGHT CALCULATION STARTED ====");
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Number of candidates: " + candidates.size());
        System.out.println("Gene tree bipartition groups: " + hashToBipartitions.size());
        
        long startTime = System.currentTimeMillis();
        
        Map<STBipartition, Double> result;
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
    private Map<STBipartition, Double> calculateWeightsSingleThread(List<STBipartition> candidates) {
        System.out.println("Starting single-threaded range-based weight calculation...");
        
        Map<STBipartition, Double> weights = new HashMap<>();
        
        // Convert candidates to range representations for efficient processing
        System.out.println("Converting candidates to range representations...");
        Map<STBipartition, RangeBipartition> candidateToRange = buildCandidateRangeMapping(candidates);
        System.out.println("Successfully mapped " + candidateToRange.size() + "/" + candidates.size() + " candidates to ranges");
        
        int processedCandidates = 0;
        
        for (STBipartition candidate : candidates) {
            double totalScore = 0.0;
            RangeBipartition candidateRange = candidateToRange.get(candidate);
            
            if (candidateRange != null) {
                // Use range-based calculation
                totalScore = calculateRangeBasedScore(candidateRange);
            } else {
                // Fallback to traditional BitSet-based calculation for unmapped candidates
                totalScore = calculateTraditionalScore(candidate);
            }
            
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
    private Map<STBipartition, Double> calculateWeightsMultiThread(List<STBipartition> candidates) {
        System.out.println("Starting multi-threaded range-based weight calculation...");
        
        Map<STBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        Threading.startThreading(numThreads);
        
        // Convert candidates to range representations
        System.out.println("Converting candidates to range representations...");
        Map<STBipartition, RangeBipartition> candidateToRange = buildCandidateRangeMapping(candidates);
        System.out.println("Successfully mapped " + candidateToRange.size() + "/" + candidates.size() + " candidates to ranges");
        
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
                        STBipartition candidate = candidates.get(j);
                        double totalScore = 0.0;
                        RangeBipartition candidateRange = candidateToRange.get(candidate);
                        
                        if (candidateRange != null) {
                            // Use range-based calculation
                            totalScore = calculateRangeBasedScore(candidateRange);
                        } else {
                            // Fallback to traditional calculation
                            totalScore = calculateTraditionalScore(candidate);
                        }
                        
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
     * Calculate score for a range bipartition against all gene tree bipartitions.
     */
    private double calculateRangeBasedScore(RangeBipartition candidateRange) {
        double totalScore = 0.0;
        
        // Iterate through gene tree bipartition groups
        for (Map.Entry<Object, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
            List<RangeBipartition> geneTreeRanges = entry.getValue();
            int frequency = geneTreeRanges.size(); // All ranges in same hash group
            
            if (!geneTreeRanges.isEmpty()) {
                RangeBipartition geneTreeRange = geneTreeRanges.get(0); // Representative
                double score = calculateRangeScore(candidateRange, geneTreeRange);
                totalScore += score * frequency;
                
                totalScoreCalculations++;
            }
        }
        
        return totalScore;
    }
    
    /**
     * Calculate score between two range bipartitions using inverse index.
     * Implements the same scoring formula as original calculateScore method.
     */
    private double calculateRangeScore(RangeBipartition range1, RangeBipartition range2) {
        // Calculate four intersection sizes: AA, AB, BA, BB
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
     * Fallback to traditional BitSet-based score calculation.
     * Used when range mapping is not available for a candidate.
     */
    private double calculateTraditionalScore(STBipartition candidate) {
        double totalScore = 0.0;
        
        // Use the existing gene tree STBipartitions from GeneTrees
        for (Map.Entry<STBipartition, Integer> entry : geneTrees.stBipartitions.entrySet()) {
            STBipartition geneTreeSTB = entry.getKey();
            int frequency = entry.getValue();
            
            double score = calculateBitSetScore(candidate, geneTreeSTB);
            totalScore += score * frequency;
            
            totalScoreCalculations++;
        }
        
        return totalScore;
    }
    
    /**
     * Traditional BitSet-based score calculation (same as original WeightCalculator).
     */
    private double calculateBitSetScore(STBipartition stb1, STBipartition stb2) {
        // First configuration: (A|B) with (X|Y)
        utils.BitSet aIntersectX = (utils.BitSet) stb1.cluster1.clone();
        aIntersectX.and(stb2.cluster1);
        int p1 = aIntersectX.cardinality();
        
        utils.BitSet bIntersectY = (utils.BitSet) stb1.cluster2.clone();
        bIntersectY.and(stb2.cluster2);
        int p2 = bIntersectY.cardinality();
        
        double score1 = 0;
        if (p1 + p2 >= 2) {
            score1 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        // Second configuration: (A|B) with (Y|X) - cross configuration
        utils.BitSet aIntersectY = (utils.BitSet) stb1.cluster1.clone();
        aIntersectY.and(stb2.cluster2);
        p1 = aIntersectY.cardinality();
        
        utils.BitSet bIntersectX = (utils.BitSet) stb1.cluster2.clone();
        bIntersectX.and(stb2.cluster1);
        p2 = bIntersectX.cardinality();
        
        double score2 = 0;
        if (p1 + p2 >= 2) {
            score2 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        totalIntersectionCalculations += 4; // 4 BitSet intersections
        
        return score1 + score2;
    }
    
    /**
     * Build mapping from STBipartitions to RangeBipartitions.
     * This is a complex operation that tries to find range representations for candidates.
     */
    private Map<STBipartition, RangeBipartition> buildCandidateRangeMapping(List<STBipartition> candidates) {
        Map<STBipartition, RangeBipartition> mapping = new HashMap<>();
        
        System.out.println("Building candidate-to-range mapping...");
        System.out.println("This may take some time for large candidate sets...");
        
        int mappedCount = 0;
        int totalCandidates = candidates.size();
        
        // For each candidate, try to find a matching range bipartition
        for (STBipartition candidate : candidates) {
            RangeBipartition matchingRange = findMatchingRange(candidate);
            if (matchingRange != null) {
                mapping.put(candidate, matchingRange);
                mappedCount++;
            }
            
            // Log progress for large datasets
            if ((mappedCount + 1) % 500 == 0 || mappedCount == totalCandidates) {
                System.out.println("Mapped " + mappedCount + "/" + totalCandidates + " candidates to ranges");
            }
        }
        
        System.out.println("Candidate mapping completed: " + mappedCount + "/" + totalCandidates + " successfully mapped");
        if (mappedCount < totalCandidates) {
            System.out.println("Unmapped candidates will use traditional BitSet calculation");
        }
        
        return mapping;
    }
    
    /**
     * Find a range bipartition that matches the given STBipartition.
     * This is an expensive operation that compares taxon sets.
     */
    private RangeBipartition findMatchingRange(STBipartition target) {
        // Convert target clusters to taxon sets for comparison
        Set<Integer> targetLeft = bitSetToTaxonSet(target.cluster1);
        Set<Integer> targetRight = bitSetToTaxonSet(target.cluster2);
        
        // Search through all range bipartitions
        for (List<RangeBipartition> ranges : hashToBipartitions.values()) {
            for (RangeBipartition range : ranges) {
                Set<Integer> rangeLeft = getRangeTaxonSet(range.geneTreeIndex, range.leftStart, range.leftEnd);
                Set<Integer> rangeRight = getRangeTaxonSet(range.geneTreeIndex, range.rightStart, range.rightEnd);
                
                // Check both orientations: (left, right) and (right, left)
                if ((targetLeft.equals(rangeLeft) && targetRight.equals(rangeRight)) ||
                    (targetLeft.equals(rangeRight) && targetRight.equals(rangeLeft))) {
                    return range;
                }
            }
        }
        
        return null; // No matching range found
    }
    
    /**
     * Convert BitSet to set of taxon IDs.
     */
    private Set<Integer> bitSetToTaxonSet(utils.BitSet bitSet) {
        Set<Integer> taxonSet = new HashSet<>();
        for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i + 1)) {
            taxonSet.add(i);
        }
        return taxonSet;
    }
    
    /**
     * Get taxon set for a range in a gene tree.
     */
    private Set<Integer> getRangeTaxonSet(int geneTreeIndex, int start, int end) {
        Set<Integer> taxonSet = new HashSet<>();
        int[][] orderings = bipartitionManager.getGeneTreeTaxaOrdering();
        
        if (geneTreeIndex >= 0 && geneTreeIndex < orderings.length) {
            int[] ordering = orderings[geneTreeIndex];
            for (int i = start; i < end && i < ordering.length; i++) {
                taxonSet.add(ordering[i]);
            }
        }
        
        return taxonSet;
    }
    
    /**
     * Memory-optimized GPU weight calculation using compact range representations.
     * This replaces the traditional BitSet-based GPU approach with O(nk) memory usage.
     */
    private Map<STBipartition, Double> calculateWeightsCompactGPU(List<STBipartition> candidates) {
        System.out.println("==== STARTING COMPACT GPU WEIGHT CALCULATION ====");
        
        try {
            // Convert candidates and gene tree bipartitions to compact representations
            System.out.println("Converting bipartitions to compact range representations...");
            
            Map<STBipartition, RangeBipartition> candidateToRange = buildCandidateRangeMapping(candidates);
            List<RangeBipartition> candidateRanges = new ArrayList<>();
            List<STBipartition> mappedCandidates = new ArrayList<>();
            
            // Build lists of successfully mapped candidates
            for (STBipartition candidate : candidates) {
                RangeBipartition range = candidateToRange.get(candidate);
                if (range != null) {
                    candidateRanges.add(range);
                    mappedCandidates.add(candidate);
                }
            }
            
            System.out.println("Successfully mapped " + candidateRanges.size() + "/" + candidates.size() + " candidates to ranges");
            
            if (candidateRanges.isEmpty()) {
                System.err.println("No candidates could be mapped to ranges - falling back to CPU calculation");
                return calculateWeightsMultiThread(candidates);
            }
            
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
            
            // Prepare GPU data structures
            int numCandidates = candidateRanges.size();
            int numGeneTreeBips = geneTreeRanges.size();
            int numTrees = inverseIndexManager.getNumTrees();
            int numTaxa = inverseIndexManager.getNumTaxa();
            
            System.out.println("==== COMPACT GPU PARAMETERS ====");
            System.out.println("Candidates: " + numCandidates);
            System.out.println("Gene tree bipartitions: " + numGeneTreeBips);
            System.out.println("Trees: " + numTrees);
            System.out.println("Taxa: " + numTaxa);
            
            // Calculate memory usage (much smaller than BitSet approach)
            long candidateMemory = (long) numCandidates * 5 * 4; // 5 ints per CompactBipartition
            long geneTreeMemory = (long) numGeneTreeBips * 5 * 4;
            long inverseIndexMemory = (long) numTrees * numTaxa * 4;
            long orderingMemory = (long) numTrees * numTaxa * 4;
            long frequencyMemory = (long) numGeneTreeBips * 4;
            long weightsMemory = (long) numCandidates * 8;
            long totalMemory = candidateMemory + geneTreeMemory + inverseIndexMemory + orderingMemory + frequencyMemory + weightsMemory;
            
            System.out.println("==== COMPACT GPU MEMORY USAGE ====");
            System.out.println("Candidate ranges: " + (candidateMemory / 1024) + " KB");
            System.out.println("Gene tree ranges: " + (geneTreeMemory / 1024) + " KB");
            System.out.println("Inverse index: " + (inverseIndexMemory / (1024 * 1024)) + " MB");
            System.out.println("Orderings: " + (orderingMemory / (1024 * 1024)) + " MB");
            System.out.println("Total GPU memory: " + (totalMemory / (1024 * 1024)) + " MB");
            System.out.println("Memory reduction vs BitSet: ~" + 
                             (((long) numCandidates + numGeneTreeBips) * numTaxa * 8 / totalMemory) + "x smaller");
            
            // Create compact bipartition arrays
            WeightCalcLib.CompactBipartition[] compactCandidates = 
                (WeightCalcLib.CompactBipartition[]) new WeightCalcLib.CompactBipartition().toArray(numCandidates);
            WeightCalcLib.CompactBipartition[] compactGeneTreeBips = 
                (WeightCalcLib.CompactBipartition[]) new WeightCalcLib.CompactBipartition().toArray(numGeneTreeBips);
            
            // Fill candidate data
            System.out.println("Preparing candidate data for GPU...");
            for (int i = 0; i < numCandidates; i++) {
                RangeBipartition range = candidateRanges.get(i);
                WeightCalcLib.CompactBipartition compact = compactCandidates[i];
                
                compact.geneTreeIndex = range.geneTreeIndex;
                compact.leftStart = range.leftStart;
                compact.leftEnd = range.leftEnd;
                compact.rightStart = range.rightStart;
                compact.rightEnd = range.rightEnd;
                compact.write();
            }
            
            // Fill gene tree data
            System.out.println("Preparing gene tree data for GPU...");
            int[] frequencyArray = new int[numGeneTreeBips];
            for (int i = 0; i < numGeneTreeBips; i++) {
                RangeBipartition range = geneTreeRanges.get(i);
                WeightCalcLib.CompactBipartition compact = compactGeneTreeBips[i];
                
                compact.geneTreeIndex = range.geneTreeIndex;
                compact.leftStart = range.leftStart;
                compact.leftEnd = range.leftEnd;
                compact.rightStart = range.rightStart;
                compact.rightEnd = range.rightEnd;
                compact.write();
                
                frequencyArray[i] = frequencies.get(i);
            }
            
            // Prepare inverse index and ordering arrays
            System.out.println("Preparing inverse index and ordering data for GPU...");
            Memory inverseIndexMem = flattenInverseIndex();
            Memory orderingMem = flattenOrderings();
            
            // Prepare results array
            double[] weights = new double[numCandidates];
            
            // Launch compact GPU kernel
            System.out.println("Launching compact GPU kernel...");
            WeightCalcLib.INSTANCE.launchCompactWeightCalculation(
                compactCandidates,
                compactGeneTreeBips,
                frequencyArray,
                weights,
                inverseIndexMem,
                orderingMem,
                numCandidates,
                numGeneTreeBips,
                numTrees,
                numTaxa
            );
            
            System.out.println("==== COMPACT GPU KERNEL COMPLETED ====");
            
            // Convert results back to Java Map
            Map<STBipartition, Double> result = new HashMap<>();
            for (int i = 0; i < numCandidates; i++) {
                result.put(mappedCandidates.get(i), weights[i]);
            }
            
            // Handle unmapped candidates with CPU fallback
            if (mappedCandidates.size() < candidates.size()) {
                System.out.println("Processing " + (candidates.size() - mappedCandidates.size()) + 
                                 " unmapped candidates with CPU fallback...");
                
                Set<STBipartition> mappedSet = new HashSet<>(mappedCandidates);
                List<STBipartition> unmappedCandidates = new ArrayList<>();
                for (STBipartition candidate : candidates) {
                    if (!mappedSet.contains(candidate)) {
                        unmappedCandidates.add(candidate);
                    }
                }
                
                // Calculate weights for unmapped candidates using traditional approach
                for (STBipartition candidate : unmappedCandidates) {
                    double weight = calculateTraditionalScore(candidate);
                    result.put(candidate, weight);
                }
            }
            
            System.out.println("Compact GPU calculation completed for " + result.size() + " candidates");
            return result;
            
        } catch (Exception e) {
            System.err.println("Compact GPU calculation failed: " + e.getMessage());
            e.printStackTrace();
            System.out.println("Falling back to CPU calculation...");
            return calculateWeightsMultiThread(candidates);
        }
    }
    
    /**
     * Flatten inverse index array for GPU transfer.
     * Format: [tree0_taxon0_pos, tree0_taxon1_pos, ..., tree1_taxon0_pos, ...]
     */
    private Memory flattenInverseIndex() {
        int[][] inverseIndex = inverseIndexManager.getInverseIndex();
        int numTrees = inverseIndex.length;
        int numTaxa = inverseIndex[0].length;
        
        Memory memory = new Memory((long) numTrees * numTaxa * 4); // 4 bytes per int
        
        for (int tree = 0; tree < numTrees; tree++) {
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                long offset = ((long) tree * numTaxa + taxon) * 4;
                memory.setInt(offset, inverseIndex[tree][taxon]);
            }
        }
        
        return memory;
    }
    
    /**
     * Flatten gene tree orderings for GPU transfer.
     * Format: [tree0_pos0_taxon, tree0_pos1_taxon, ..., tree1_pos0_taxon, ...]
     */
    private Memory flattenOrderings() {
        int[][] orderings = inverseIndexManager.getGeneTreeOrderings();
        int numTrees = orderings.length;
        int numTaxa = orderings[0].length;
        
        Memory memory = new Memory((long) numTrees * numTaxa * 4); // 4 bytes per int
        
        for (int tree = 0; tree < numTrees; tree++) {
            for (int pos = 0; pos < numTaxa; pos++) {
                long offset = ((long) tree * numTaxa + pos) * 4;
                memory.setInt(offset, orderings[tree][pos]);
            }
        }
        
        return memory;
    }
    
    /**
     * Get processing statistics for performance analysis.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Memory-Optimized Weight Calculator Statistics:\n");
        sb.append("  Total processing time: ").append(totalProcessingTime).append(" ms\n");
        sb.append("  Total score calculations: ").append(totalScoreCalculations).append("\n");
        sb.append("  Total intersection calculations: ").append(totalIntersectionCalculations).append("\n");
        sb.append("  Range bipartition groups: ").append(hashToBipartitions.size()).append("\n");
        
        if (totalScoreCalculations > 0) {
            sb.append("  Average intersections per score: ").append(totalIntersectionCalculations / (double) totalScoreCalculations).append("\n");
            sb.append("  Average time per score: ").append(totalProcessingTime / (double) totalScoreCalculations).append(" ms\n");
        }
        
        return sb.toString();
    }
}
