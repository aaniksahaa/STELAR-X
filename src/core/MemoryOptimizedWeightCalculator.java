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
    private final RangeBipartition[] geneTreeRanges;
    private final int[] geneTreeRangeFrequencies;
    private final Map<RangeBipartition, Integer> rangeIndex;
    private final int bitsetWords;
    private final long[] leftRangeBits;
    private final long[] rightRangeBits;

    private static final long CPU_BITSET_MEMORY_CAP = 64L * 1024 * 1024;
    
    // Statistics for performance monitoring
    private long totalScoreCalculations = 0;
    private long totalIntersectionCalculations = 0;
    private long totalProcessingTime = 0;
    
    public MemoryOptimizedWeightCalculator(GeneTrees geneTrees) {
        System.out.println("==== INITIALIZING MEMORY-OPTIMIZED WEIGHT CALCULATOR ====");
        
        this.geneTrees = geneTrees;
        
        System.out.println("Reusing preprocessed range bipartitions...");
        this.bipartitionManager = geneTrees.getBipartitionManager();
        
        System.out.println("Creating InverseIndexManager...");
        this.inverseIndexManager = new InverseIndexManager(
            bipartitionManager.getGeneTreeTaxaOrdering(), geneTrees.realTaxaCount);
        
        this.hashToBipartitions = bipartitionManager.getHashToBipartitions();
        this.geneTreeRanges = new RangeBipartition[geneTrees.rangeBipartitions.size()];
        this.geneTreeRangeFrequencies = new int[geneTreeRanges.length];
        this.rangeIndex = new HashMap<>(hashMapCapacity(geneTreeRanges.length));
        int rangeIndex = 0;
        for (Map.Entry<RangeBipartition, Integer> entry : geneTrees.rangeBipartitions.entrySet()) {
            geneTreeRanges[rangeIndex] = entry.getKey();
            geneTreeRangeFrequencies[rangeIndex] = entry.getValue();
            this.rangeIndex.put(entry.getKey(), rangeIndex);
            rangeIndex++;
        }

        int words = (geneTrees.realTaxaCount + Long.SIZE - 1) / Long.SIZE;
        long bitsetBytes = 2L * geneTreeRanges.length * words * Long.BYTES;
        double averageUnionSize = averageUnionSize(geneTreeRanges);
        boolean useBitsets = Config.COMPUTATION_MODE != Config.ComputationMode.GPU_PARALLEL
            && words > 0
            && bitsetBytes <= CPU_BITSET_MEMORY_CAP
            && 4.0 * words <= averageUnionSize;
        if (useBitsets) {
            this.bitsetWords = words;
            this.leftRangeBits = new long[geneTreeRanges.length * words];
            this.rightRangeBits = new long[geneTreeRanges.length * words];
            buildRangeBitsets();
            System.out.println("CPU intersection strategy: bitset popcount ("
                + (bitsetBytes / 1024) + " KB resident)");
        } else {
            this.bitsetWords = 0;
            this.leftRangeBits = null;
            this.rightRangeBits = null;
            System.out.println("CPU intersection strategy: fused inverse-index range walk");
        }
        
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
        
        double[] scores = new double[candidates.size()];
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
                        scores[j] = calculateRangeBasedScore(candidate);
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
        
        // Each worker writes to a disjoint score-array slice. Populate the result map
        // once, serially, instead of contending on ConcurrentHashMap in the hot loop.
        Map<RangeBipartition, Double> weights = new HashMap<>(hashMapCapacity(candidates.size()));
        for (int i = 0; i < candidates.size(); i++) {
            weights.put(candidates.get(i), scores[i]);
        }

        System.out.println("Multi-threaded calculation completed");
        return weights;
    }

    private static int hashMapCapacity(int expectedSize) {
        if (expectedSize < 3) return expectedSize + 1;
        return expectedSize < (1 << 29) ? (int) (expectedSize / 0.75f) + 1 : Integer.MAX_VALUE;
    }

    private static double averageUnionSize(RangeBipartition[] ranges) {
        if (ranges.length == 0) return 0.0;
        long total = 0;
        for (RangeBipartition range : ranges) {
            total += range.leftSize() + range.rightSize();
        }
        return total / (double) ranges.length;
    }

    private void buildRangeBitsets() {
        int[][] orderings = bipartitionManager.getGeneTreeTaxaOrdering();
        for (int i = 0; i < geneTreeRanges.length; i++) {
            RangeBipartition range = geneTreeRanges[i];
            int[] ordering = orderings[range.geneTreeIndex];
            int base = i * bitsetWords;
            setRangeBits(leftRangeBits, base, ordering, range.leftStart, range.leftEnd);
            setRangeBits(rightRangeBits, base, ordering, range.rightStart, range.rightEnd);
        }
    }

    private static void setRangeBits(long[] bits, int base, int[] ordering, int start, int end) {
        for (int pos = start; pos < end; pos++) {
            int taxon = ordering[pos];
            bits[base + (taxon >>> 6)] |= 1L << (taxon & 63);
        }
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
        if (bitsetWords != 0) {
            Integer candidateIndex = rangeIndex.get(candidate);
            if (candidateIndex != null) {
                return calculateBitsetScore(candidateIndex);
            }
        }

        double totalScore = 0.0;
        int[] intersections = new int[4];

        // Dense arrays avoid Map.Entry iteration in this quadratic hot loop. The
        // fused call fills AA, AB, BA, BB in one traversal of the smaller union.
        for (int i = 0; i < geneTreeRanges.length; i++) {
            RangeBipartition geneTreeRange = geneTreeRanges[i];
            inverseIndexManager.getBipartitionIntersections(
                candidate.geneTreeIndex,
                candidate.leftStart, candidate.leftEnd,
                candidate.rightStart, candidate.rightEnd,
                geneTreeRange.geneTreeIndex,
                geneTreeRange.leftStart, geneTreeRange.leftEnd,
                geneTreeRange.rightStart, geneTreeRange.rightEnd,
                intersections);

            int aa = intersections[0];
            int ab = intersections[1];
            int ba = intersections[2];
            int bb = intersections[3];
            double score = 0.0;
            int aligned = aa + bb;
            if (aligned >= 2) score += aa * (double) bb * (aligned - 2) / 2.0;
            int crossed = ab + ba;
            if (crossed >= 2) score += ab * (double) ba * (crossed - 2) / 2.0;

            totalScore += score * geneTreeRangeFrequencies[i];
        }
        totalScoreCalculations += geneTreeRanges.length;
        totalIntersectionCalculations += 4L * geneTreeRanges.length;
        
        return totalScore;
    }

    private double calculateBitsetScore(int candidateIndex) {
        double totalScore = 0.0;
        int candidateBase = candidateIndex * bitsetWords;

        for (int i = 0; i < geneTreeRanges.length; i++) {
            int geneBase = i * bitsetWords;
            int aa = 0, ab = 0, ba = 0, bb = 0;
            for (int word = 0; word < bitsetWords; word++) {
                long candidateLeft = leftRangeBits[candidateBase + word];
                long candidateRight = rightRangeBits[candidateBase + word];
                long geneLeft = leftRangeBits[geneBase + word];
                long geneRight = rightRangeBits[geneBase + word];
                aa += Long.bitCount(candidateLeft & geneLeft);
                ab += Long.bitCount(candidateLeft & geneRight);
                ba += Long.bitCount(candidateRight & geneLeft);
                bb += Long.bitCount(candidateRight & geneRight);
            }

            double score = 0.0;
            int aligned = aa + bb;
            if (aligned >= 2) score += aa * (double) bb * (aligned - 2) / 2.0;
            int crossed = ab + ba;
            if (crossed >= 2) score += ab * (double) ba * (crossed - 2) / 2.0;
            totalScore += score * geneTreeRangeFrequencies[i];
        }
        totalScoreCalculations += geneTreeRanges.length;
        totalIntersectionCalculations += 4L * geneTreeRanges.length;
        return totalScore;
    }
    
    // ========================================================================
    // MIXED BIPARTITION SUPPORT
    // For cross-tree recombination where left and right sides may come from
    // different gene trees.
    // ========================================================================
    
    /**
     * Calculate score between a MixedBipartition and a RangeBipartition (gene tree).
     * 
     * Same intersection logic as calculateRangeScore, but uses separate tree indices
     * for left and right sides of the mixed bipartition.
     * 
     * @param mixed The mixed bipartition (sides may be from different trees)
     * @param geneTree The gene tree bipartition
     * @return The score contribution
     */
    private double calculateMixedScore(MixedBipartition mixed, RangeBipartition geneTree) {
        // Calculate four intersection sizes: AA, AB, BA, BB
        // Key difference: use leftTreeIndex for left side, rightTreeIndex for right side
        int aa = inverseIndexManager.getRangeIntersectionSize(
            mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
            geneTree.geneTreeIndex, geneTree.leftStart, geneTree.leftEnd);
            
        int bb = inverseIndexManager.getRangeIntersectionSize(
            mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
            geneTree.geneTreeIndex, geneTree.rightStart, geneTree.rightEnd);
            
        int ab = inverseIndexManager.getRangeIntersectionSize(
            mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
            geneTree.geneTreeIndex, geneTree.rightStart, geneTree.rightEnd);
            
        int ba = inverseIndexManager.getRangeIntersectionSize(
            mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
            geneTree.geneTreeIndex, geneTree.leftStart, geneTree.leftEnd);
        
        totalIntersectionCalculations += 4;
        
        // Same scoring formula as RangeBipartition
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
     * Calculate total weight for a single MixedBipartition.
     * Sums the score contribution from all gene tree bipartitions.
     * 
     * @param mixed The mixed bipartition
     * @return The total weight
     */
    public double calculateMixedWeight(MixedBipartition mixed) {
        double totalScore = 0.0;
        
        for (Map.Entry<RangeBipartition, Integer> entry : geneTrees.rangeBipartitions.entrySet()) {
            RangeBipartition geneTreeRange = entry.getKey();
            int frequency = entry.getValue();
            
            double score = calculateMixedScore(mixed, geneTreeRange);
            totalScore += score * frequency;
            
            totalScoreCalculations++;
        }
        
        return totalScore;
    }
    
    /**
     * Calculate weights for a list of MixedBipartitions.
     * Uses the configured computation mode (CPU single/parallel/GPU).
     * 
     * @param mixedCandidates List of mixed bipartitions
     * @return Map of mixed bipartition to weight
     */
    public Map<MixedBipartition, Double> calculateMixedWeights(List<MixedBipartition> mixedCandidates) {
        System.out.println("==== MIXED BIPARTITION WEIGHT CALCULATION STARTED ====");
        System.out.println("Number of mixed candidates: " + mixedCandidates.size());
        
        long startTime = System.currentTimeMillis();
        
        Map<MixedBipartition, Double> result;
        switch (Config.COMPUTATION_MODE) {
            case CPU_SINGLE:
                result = calculateMixedWeightsSingleThread(mixedCandidates);
                break;
            case CPU_PARALLEL:
                result = calculateMixedWeightsMultiThread(mixedCandidates);
                break;
            case GPU_PARALLEL:
                System.out.println("Using GPU for mixed bipartition weight calculation");
                result = calculateMixedWeightsGPU(mixedCandidates);
                break;
            default:
                throw new IllegalStateException("Unknown computation mode: " + Config.COMPUTATION_MODE);
        }
        
        long endTime = System.currentTimeMillis();
        
        System.out.println("==== MIXED BIPARTITION WEIGHT CALCULATION COMPLETED ====");
        System.out.println("Processing time: " + (endTime - startTime) + " ms");
        System.out.println("Mixed bipartitions processed: " + result.size());
        
        return result;
    }
    
    /**
     * Single-threaded weight calculation for mixed bipartitions.
     */
    private Map<MixedBipartition, Double> calculateMixedWeightsSingleThread(List<MixedBipartition> candidates) {
        Map<MixedBipartition, Double> weights = new HashMap<>();
        
        int processed = 0;
        for (MixedBipartition mixed : candidates) {
            double weight = calculateMixedWeight(mixed);
            weights.put(mixed, weight);
            
            processed++;
            if (processed % 1000 == 0 || processed == candidates.size()) {
                System.out.println("Processed " + processed + "/" + candidates.size() + " mixed bipartitions");
            }
        }
        
        return weights;
    }
    
    /**
     * Multi-threaded weight calculation for mixed bipartitions.
     */
    private Map<MixedBipartition, Double> calculateMixedWeightsMultiThread(List<MixedBipartition> candidates) {
        Map<MixedBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        Threading.startThreading(numThreads);
        
        int chunkSize = Math.max(1, (candidates.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (candidates.size() + chunkSize - 1) / chunkSize);
        CountDownLatch latch = new CountDownLatch(actualThreads);
        
        System.out.println("Using " + actualThreads + " threads for mixed bipartition weight calculation");
        
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, candidates.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    if (startIdx >= endIdx || startIdx >= candidates.size()) {
                        return;
                    }
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        MixedBipartition mixed = candidates.get(j);
                        double weight = calculateMixedWeight(mixed);
                        weights.put(mixed, weight);
                    }
                    
                    System.out.println("Thread " + threadId + " completed " + (endIdx - startIdx) + " mixed bipartitions");
                    
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Mixed weight calculation interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        return weights;
    }
    
    /**
     * GPU weight calculation for mixed bipartitions using MixedCompactBipartition kernel.
     */
    private Map<MixedBipartition, Double> calculateMixedWeightsGPU(List<MixedBipartition> candidates) {
        System.out.println("==== STARTING MIXED BIPARTITION GPU CALCULATION ====");
        
        try {
            // Check if CUDA library is available
            try {
                WeightCalcLib.INSTANCE.toString();
                System.out.println("CUDA library loaded successfully for mixed bipartitions");
            } catch (UnsatisfiedLinkError e) {
                System.err.println("CUDA library not found: " + e.getMessage());
                System.out.println("Falling back to CPU calculation for mixed bipartitions...");
                return calculateMixedWeightsMultiThread(candidates);
            }
            
            int numCandidates = candidates.size();
            int numTrees = inverseIndexManager.getNumTrees();
            int numTaxa = inverseIndexManager.getNumTaxa();
            
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
            
            int numGeneTreeBips = geneTreeRanges.size();
            
            System.out.println("Mixed GPU Parameters:");
            System.out.println("  Mixed candidates: " + numCandidates);
            System.out.println("  Gene tree bipartitions: " + numGeneTreeBips);
            System.out.println("  Trees: " + numTrees);
            System.out.println("  Taxa: " + numTaxa);
            
            // Convert mixed candidates to GPU format using contiguous memory allocation
            WeightCalcLib.MixedCompactBipartition mixedTemplate = new WeightCalcLib.MixedCompactBipartition();
            WeightCalcLib.MixedCompactBipartition[] mixedArray = 
                (WeightCalcLib.MixedCompactBipartition[]) mixedTemplate.toArray(numCandidates);
            
            for (int i = 0; i < numCandidates; i++) {
                MixedBipartition mixed = candidates.get(i);
                mixedArray[i].leftTreeIndex = mixed.leftTreeIndex;
                mixedArray[i].leftStart = mixed.leftStart;
                mixedArray[i].leftEnd = mixed.leftEnd;
                mixedArray[i].rightTreeIndex = mixed.rightTreeIndex;
                mixedArray[i].rightStart = mixed.rightStart;
                mixedArray[i].rightEnd = mixed.rightEnd;
            }
            
            // Convert gene tree ranges to GPU compact format
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
            
            // Launch mixed bipartition GPU kernel
            System.out.println("Launching mixed bipartition GPU kernel...");
            System.out.flush(); // Flush to ensure we see logs before native call
            
            try {
                WeightCalcLib.INSTANCE.launchMixedWeightCalculation(
                    mixedArray,
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
            } catch (UnsatisfiedLinkError e) {
                System.err.println("ERROR: launchMixedWeightCalculation not found in CUDA library!");
                System.err.println("The CUDA library needs to be rebuilt with the new mixed bipartition kernel.");
                System.err.println("Run: ./build.sh to rebuild the CUDA library");
                System.err.println("JNA Error: " + e.getMessage());
                throw e;
            } catch (Exception e) {
                System.err.println("ERROR calling launchMixedWeightCalculation: " + e.getMessage());
                e.printStackTrace();
                throw e;
            }
            
            System.out.println("==== MIXED BIPARTITION GPU KERNEL COMPLETED ====");
            
            // Convert results back to Java Map
            Map<MixedBipartition, Double> result = new HashMap<>();
            for (int i = 0; i < numCandidates; i++) {
                result.put(candidates.get(i), weights[i]);
            }
            
            return result;
            
        } catch (Exception e) {
            System.err.println("Mixed bipartition GPU calculation failed: " + e.getMessage());
            e.printStackTrace();
            System.out.println("Falling back to CPU calculation...");
            return calculateMixedWeightsMultiThread(candidates);
        }
    }
    
    /**
     * Get the inverse index manager (for external use).
     */
    public InverseIndexManager getInverseIndexManager() {
        return inverseIndexManager;
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
