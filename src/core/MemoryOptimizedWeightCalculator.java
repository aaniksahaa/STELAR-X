package core;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;

import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
import com.sun.jna.Structure;

import preprocessing.MemoryOptimizedGeneTrees;
import tree.RangeSTBipartition;
import utils.Config;
import utils.Threading;

/**
 * Memory-optimized weight calculator using range-based bipartitions.
 * Instead of bitset operations, uses array-based intersection counting.
 */
public class MemoryOptimizedWeightCalculator {
    
    private MemoryOptimizedGeneTrees geneTrees;
    
    // JNA interface for CUDA library with range-based bipartitions
    public interface RangeWeightCalcLib extends Library {
        RangeWeightCalcLib INSTANCE = Native.load("range_weight_calc", RangeWeightCalcLib.class);
        
        // Structure for range-based bipartition
        @Structure.FieldOrder({"geneTreeIndex", "startPos", "endPos"})
        public static class RangeBipartition extends Structure {
            public int geneTreeIndex;
            public int startPos;
            public int endPos;
            
            public RangeBipartition() {
                super();
            }
            
            public RangeBipartition(Pointer p) {
                super(p);
            }
        }
        
        void launchRangeWeightCalculation(
            RangeBipartition[] candidates,
            RangeBipartition[] geneTreeBips,
            int[] frequencies,
            double[] weights,
            int[][] geneTreeOrderings,
            int[] orderingLengths,
            int numCandidates,
            int numGeneTreeBips,
            int numGeneTrees,
            int maxTaxaCount
        );
    }
    
    public MemoryOptimizedWeightCalculator(MemoryOptimizedGeneTrees geneTrees) {
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
                return calculateWeightsGPU(candidates);
            default:
                throw new IllegalStateException("Unknown computation mode: " + Config.COMPUTATION_MODE);
        }
    }
    
    private Map<RangeSTBipartition, Double> calculateWeightsSingleThread(List<RangeSTBipartition> candidates) {
        Map<RangeSTBipartition, Double> weights = new HashMap<>();
        int[][] orderings = geneTrees.getGeneTreeOrderings();
        
        for (RangeSTBipartition candidate : candidates) {
            double totalScore = 0.0;
            
            for (Map.Entry<RangeSTBipartition, Integer> entry : geneTrees.rangeSTBipartitions.entrySet()) {
                RangeSTBipartition geneTreeSTB = entry.getKey();
                int frequency = entry.getValue();
                
                double score = calculateScore(candidate, geneTreeSTB, orderings);
                totalScore += score * frequency;
            }
            
            weights.put(candidate, totalScore);
        }
        
        return weights;
    }
    
    private Map<RangeSTBipartition, Double> calculateWeightsMultiThread(List<RangeSTBipartition> candidates) {
        Map<RangeSTBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        int[][] orderings = geneTrees.getGeneTreeOrderings();
        
        Threading.startThreading(numThreads);
        
        int chunkSize = (candidates.size() + numThreads - 1) / numThreads;
        CountDownLatch latch = new CountDownLatch(numThreads);
        
        for (int i = 0; i < numThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, candidates.size());
            
            Threading.execute(() -> {
                try {
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
        
        return weights;
    }
    
    private Map<RangeSTBipartition, Double> calculateWeightsGPU(List<RangeSTBipartition> candidates) {
        try {
            System.out.println("Starting GPU range weight calculation...");
            
            int numCandidates = candidates.size();
            int numGeneTreeBips = geneTrees.rangeSTBipartitions.size();
            int numGeneTrees = geneTrees.geneTrees.size();
            int maxTaxaCount = geneTrees.realTaxaCount;
            
            System.out.println("Number of candidates: " + numCandidates);
            System.out.println("Number of gene tree bipartitions: " + numGeneTreeBips);
            System.out.println("Number of gene trees: " + numGeneTrees);
            System.out.println("Max taxa count: " + maxTaxaCount);
            
            // Allocate arrays for CUDA
            RangeWeightCalcLib.RangeBipartition[] cudaCandidates = 
                (RangeWeightCalcLib.RangeBipartition[]) new RangeWeightCalcLib.RangeBipartition().toArray(numCandidates);
            RangeWeightCalcLib.RangeBipartition[] cudaGeneTreeBips = 
                (RangeWeightCalcLib.RangeBipartition[]) new RangeWeightCalcLib.RangeBipartition().toArray(numGeneTreeBips);
            int[] frequencies = new int[numGeneTreeBips];
            double[] weights = new double[numCandidates];
            
            // Prepare gene tree orderings for GPU
            int[][] orderings = geneTrees.getGeneTreeOrderings();
            int[] orderingLengths = new int[numGeneTrees];
            for (int i = 0; i < numGeneTrees; i++) {
                orderingLengths[i] = orderings[i].length;
            }
            
            // Convert candidates
            for (int i = 0; i < numCandidates; i++) {
                RangeSTBipartition rangeBip = candidates.get(i);
                RangeWeightCalcLib.RangeBipartition cudaBip = cudaCandidates[i];
                
                cudaBip.geneTreeIndex = rangeBip.geneTreeIndex;
                cudaBip.startPos = rangeBip.startPos;
                cudaBip.endPos = rangeBip.endPos;
                cudaBip.write();
            }
            
            // Convert gene tree bipartitions
            int idx = 0;
            for (Map.Entry<RangeSTBipartition, Integer> entry : geneTrees.rangeSTBipartitions.entrySet()) {
                RangeSTBipartition rangeBip = entry.getKey();
                RangeWeightCalcLib.RangeBipartition cudaBip = cudaGeneTreeBips[idx];
                
                cudaBip.geneTreeIndex = rangeBip.geneTreeIndex;
                cudaBip.startPos = rangeBip.startPos;
                cudaBip.endPos = rangeBip.endPos;
                cudaBip.write();
                
                frequencies[idx] = entry.getValue();
                idx++;
            }
            
            System.out.println("Launching CUDA kernel for range bipartitions...");
            RangeWeightCalcLib.INSTANCE.launchRangeWeightCalculation(
                cudaCandidates,
                cudaGeneTreeBips,
                frequencies,
                weights,
                orderings,
                orderingLengths,
                numCandidates,
                numGeneTreeBips,
                numGeneTrees,
                maxTaxaCount
            );
            
            System.out.println("CUDA kernel completed");
            
            // Convert results back to Java Map
            Map<RangeSTBipartition, Double> result = new HashMap<>();
            for (int i = 0; i < numCandidates; i++) {
                result.put(candidates.get(i), weights[i]);
            }
            
            // Print sample weights
            System.out.println("Sample range-based weights:");
            int count = 0;
            for (Map.Entry<RangeSTBipartition, Double> entry : result.entrySet()) {
                if (count < 5) {
                    System.out.println("Weight " + count + ": " + entry.getValue() + 
                                     " for " + entry.getKey());
                    count++;
                } else {
                    break;
                }
            }
            
            return result;
            
        } catch (Exception e) {
            System.err.println("GPU computation failed, falling back to CPU: " + e.getMessage());
            e.printStackTrace();
            return calculateWeightsMultiThread(candidates);
        }
    }
}
