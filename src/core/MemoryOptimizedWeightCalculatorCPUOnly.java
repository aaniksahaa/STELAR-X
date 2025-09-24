package core;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;

import preprocessing.MemoryOptimizedGeneTrees;
import tree.RangeSTBipartition;
import utils.Config;
import utils.Threading;

/**
 * CPU-only memory-optimized weight calculator using range-based bipartitions.
 * This version doesn't use GPU/JNA dependencies and focuses on CPU computation only.
 */
public class MemoryOptimizedWeightCalculatorCPUOnly {
    
    private MemoryOptimizedGeneTrees geneTrees;
    
    public MemoryOptimizedWeightCalculatorCPUOnly(MemoryOptimizedGeneTrees geneTrees) {
        this.geneTrees = geneTrees;
    }
    
    /**
     * Calculates score between two range bipartitions using the same logic as BitSet version.
     */
    private double calculateScore(RangeSTBipartition stb1, RangeSTBipartition stb2, 
                                 int[][] geneTreeOrderings) {
        // Get taxa arrays for both bipartitions and their complements
        int[] stb1_cluster1 = stb1.getTaxaIds(geneTreeOrderings);
        int[] stb1_cluster2 = stb1.getComplementTaxaIds(geneTreeOrderings);
        int[] stb2_cluster1 = stb2.getTaxaIds(geneTreeOrderings);
        int[] stb2_cluster2 = stb2.getComplementTaxaIds(geneTreeOrderings);
        
        // First configuration: (A|B) with (X|Y)
        int p1 = calculateIntersectionSizeArrays(stb1_cluster1, stb2_cluster1);
        int p2 = calculateIntersectionSizeArrays(stb1_cluster2, stb2_cluster2);
        
        double score1 = 0;
        if (p1 + p2 >= 2) {
            score1 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        // Second configuration: (A|B) with (Y|X) - cross configuration
        p1 = calculateIntersectionSizeArrays(stb1_cluster1, stb2_cluster2);
        p2 = calculateIntersectionSizeArrays(stb1_cluster2, stb2_cluster1);
        
        double score2 = 0;
        if (p1 + p2 >= 2) {
            score2 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        return score1 + score2;
    }
    
    /**
     * Helper method to calculate intersection size between two sorted arrays.
     */
    private int calculateIntersectionSizeArrays(int[] arr1, int[] arr2) {
        // Sort arrays if not already sorted
        java.util.Arrays.sort(arr1);
        java.util.Arrays.sort(arr2);
        
        int count = 0;
        int i = 0, j = 0;
        
        while (i < arr1.length && j < arr2.length) {
            if (arr1[i] == arr2[j]) {
                count++;
                i++;
                j++;
            } else if (arr1[i] < arr2[j]) {
                i++;
            } else {
                j++;
            }
        }
        
        return count;
    }
    
    public Map<RangeSTBipartition, Double> calculateWeights(List<RangeSTBipartition> candidates) {
        switch (Config.COMPUTATION_MODE) {
            case CPU_SINGLE:
                return calculateWeightsSingleThread(candidates);
            case CPU_PARALLEL:
                return calculateWeightsMultiThread(candidates);
            case GPU_PARALLEL:
                System.out.println("GPU mode not supported in CPU-only version, falling back to CPU_PARALLEL");
                return calculateWeightsMultiThread(candidates);
            default:
                throw new IllegalStateException("Unknown computation mode: " + Config.COMPUTATION_MODE);
        }
    }
    
    private Map<RangeSTBipartition, Double> calculateWeightsSingleThread(List<RangeSTBipartition> candidates) {
        Map<RangeSTBipartition, Double> weights = new HashMap<>();
        int[][] orderings = geneTrees.getGeneTreeOrderings();
        
        System.out.println("Calculating weights for " + candidates.size() + " candidates (single-threaded)...");
        
        for (int i = 0; i < candidates.size(); i++) {
            RangeSTBipartition candidate = candidates.get(i);
            double totalScore = 0.0;
            
            for (Map.Entry<RangeSTBipartition, Integer> entry : geneTrees.rangeSTBipartitions.entrySet()) {
                RangeSTBipartition geneTreeSTB = entry.getKey();
                int frequency = entry.getValue();
                
                double score = calculateScore(candidate, geneTreeSTB, orderings);
                totalScore += score * frequency;
            }
            
            weights.put(candidate, totalScore);
            
            // Progress reporting
            if ((i + 1) % 10 == 0 || (i + 1) == candidates.size()) {
                System.out.printf("Processed %d/%d candidates (%.1f%%)\n", 
                    i + 1, candidates.size(), 100.0 * (i + 1) / candidates.size());
            }
        }
        
        return weights;
    }
    
    private Map<RangeSTBipartition, Double> calculateWeightsMultiThread(List<RangeSTBipartition> candidates) {
        Map<RangeSTBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        int[][] orderings = geneTrees.getGeneTreeOrderings();
        
        System.out.println("Calculating weights for " + candidates.size() + 
                          " candidates using " + numThreads + " threads...");
        
        Threading.startThreading(numThreads);
        
        int chunkSize = (candidates.size() + numThreads - 1) / numThreads;
        CountDownLatch latch = new CountDownLatch(numThreads);
        
        for (int i = 0; i < numThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, candidates.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    System.out.println("Thread " + threadId + " processing candidates " + 
                                     startIdx + " to " + (endIdx - 1));
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        RangeSTBipartition candidate = candidates.get(j);
                        double totalScore = 0.0;
                        
                        // Create local copy to avoid concurrent modification
                        Map<RangeSTBipartition, Integer> localBipartitions = 
                            new HashMap<>(geneTrees.rangeSTBipartitions);
                        
                        for (Map.Entry<RangeSTBipartition, Integer> entry : localBipartitions.entrySet()) {
                            RangeSTBipartition geneTreeSTB = entry.getKey();
                            int frequency = entry.getValue();
                            
                            double score = calculateScore(candidate, geneTreeSTB, orderings);
                            totalScore += score * frequency;
                        }
                        
                        weights.put(candidate, totalScore);
                    }
                    
                    System.out.println("Thread " + threadId + " completed processing " + 
                                     (endIdx - startIdx) + " candidates");
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
        
        System.out.println("Weight calculation completed for all " + candidates.size() + " candidates");
        return weights;
    }
}
