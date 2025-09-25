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
import com.sun.jna.Memory;

import preprocessing.GeneTrees;
import tree.STBipartition;
import utils.BitSet;
import utils.Config;
import utils.Threading;

// The parallelism is mainly happening here
public class WeightCalculator {
    
    private GeneTrees geneTrees;
    
    // JNA interface for CUDA library
    public interface WeightCalcLib extends Library {
        WeightCalcLib INSTANCE = Native.load("weight_calc", WeightCalcLib.class);
        
        // Structure to match CUDA Bipartition struct
        @Structure.FieldOrder({"cluster1", "cluster2", "bitsetSize"})
        public static class Bipartition extends Structure {
            public Pointer cluster1;
            public Pointer cluster2;
            public int bitsetSize;
            
            public Bipartition() {
                super();
            }
            
            public Bipartition(Pointer p) {
                super(p);
            }
        }
        
        void launchWeightCalculation(
            Bipartition[] candidates,
            Bipartition[] geneTreeBips,
            int[] frequencies,
            double[] weights,
            int numCandidates,
            int numGeneTreeBips,
            int bitsetSize
        );
    }
    
    public WeightCalculator(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
    }
    
    private double calculateScore(STBipartition stb1, STBipartition stb2) {
        // First configuration: (A|B) with (X|Y)
        BitSet aIntersectX = (BitSet) stb1.cluster1.clone();
        aIntersectX.and(stb2.cluster1);
        int p1 = aIntersectX.cardinality();
        
        BitSet bIntersectY = (BitSet) stb1.cluster2.clone();
        bIntersectY.and(stb2.cluster2);
        int p2 = bIntersectY.cardinality();
        
        double score1 = 0;
        if (p1 + p2 >= 2) {
            score1 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        // Second configuration: (A|B) with (Y|X) - cross configuration
        BitSet aIntersectY = (BitSet) stb1.cluster1.clone();
        aIntersectY.and(stb2.cluster2);
        p1 = aIntersectY.cardinality();
        
        BitSet bIntersectX = (BitSet) stb1.cluster2.clone();
        bIntersectX.and(stb2.cluster1);
        p2 = bIntersectX.cardinality();
        
        double score2 = 0;
        if (p1 + p2 >= 2) {
            score2 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        return score1 + score2;
    }
    
    public Map<STBipartition, Double> calculateWeights(List<STBipartition> candidates) {
        System.out.println("==== WEIGHT CALCULATION STARTED ====");
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Number of candidates: " + (candidates != null ? candidates.size() : "NULL"));
        System.out.println("Gene trees available: " + (geneTrees != null ? "YES" : "NO"));
        if (geneTrees != null) {
            System.out.println("Gene trees stBipartitions count: " + 
                             (geneTrees.stBipartitions != null ? geneTrees.stBipartitions.size() : "NULL"));
            System.out.println("Gene trees realTaxaCount: " + geneTrees.realTaxaCount);
        }
        
        switch (Config.COMPUTATION_MODE) {
            case CPU_SINGLE:
                System.out.println("Using CPU_SINGLE mode");
                return calculateWeightsSingleThread(candidates);
            case CPU_PARALLEL:
                System.out.println("Using CPU_PARALLEL mode");
                return calculateWeightsMultiThread(candidates);
            case GPU_PARALLEL:
                System.out.println("Using GPU_PARALLEL mode");
                return calculateWeightsGPU(candidates);
            default:
                System.out.println("Using default mode (throwing exception)");
                throw new IllegalStateException("Unknown computation mode: " + Config.COMPUTATION_MODE);
        }
    }
    
    private Map<STBipartition, Double> calculateWeightsSingleThread(List<STBipartition> candidates) {
        Map<STBipartition, Double> weights = new HashMap<>();
        
        for (STBipartition candidate : candidates) {
            double totalScore = 0.0;
            
            for (Map.Entry<STBipartition, Integer> entry : geneTrees.stBipartitions.entrySet()) {
                STBipartition geneTreeSTB = entry.getKey();
                int frequency = entry.getValue();
                
                double score = calculateScore(candidate, geneTreeSTB);
                totalScore += score * frequency;
            }
            
            weights.put(candidate, totalScore);
        }
        
        return weights;
    }
    
    private Map<STBipartition, Double> calculateWeightsMultiThread(List<STBipartition> candidates) {
        Map<STBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        Threading.startThreading(numThreads);
        
        // Calculate optimal number of threads to avoid invalid ranges
        int chunkSize = Math.max(1, (candidates.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (candidates.size() + chunkSize - 1) / chunkSize);
        CountDownLatch latch = new CountDownLatch(actualThreads);
        
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, candidates.size());
            
            Threading.execute(() -> {
                try {
                    // Validate range before processing
                    if (startIdx >= endIdx || startIdx >= candidates.size()) {
                        return;
                    }
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        STBipartition candidate = candidates.get(j);
                        double totalScore = 0.0;
                        
                        Map<STBipartition, Integer> localBipartitions = new HashMap<>(geneTrees.stBipartitions);
                        
                        for (Map.Entry<STBipartition, Integer> entry : localBipartitions.entrySet()) {
                            STBipartition geneTreeSTB = entry.getKey();
                            int frequency = entry.getValue();
                            
                            double score = calculateScore(candidate, geneTreeSTB);
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
    
    private Map<STBipartition, Double> calculateWeightsGPU(List<STBipartition> candidates) {
        try {
            System.out.println("==== STARTING GPU WEIGHT CALCULATION ====");
            
            // Convert Java data structures to CUDA-compatible format
            int numCandidates = candidates.size();
            int numGeneTreeBips = geneTrees.stBipartitions.size();
            int bitsetSize = (geneTrees.realTaxaCount + 63) / 64; // Size in 64-bit words
            
            System.out.println("==== CUDA ALLOCATION PARAMETERS ====");
            System.out.println("numCandidates: " + numCandidates);
            System.out.println("numGeneTreeBips: " + numGeneTreeBips);
            System.out.println("bitsetSize (in 64-bit words): " + bitsetSize);
            System.out.println("geneTrees.realTaxaCount: " + geneTrees.realTaxaCount);
            
            // Log potential allocation sizes
            long candidateMemorySize = (long) numCandidates * bitsetSize * 8L;
            long geneTreeMemorySize = (long) numGeneTreeBips * bitsetSize * 8L;
            System.out.println("==== MEMORY ALLOCATION SIZES ====");
            System.out.println("candidateCluster1 Memory: " + candidateMemorySize + " bytes");
            System.out.println("candidateCluster2 Memory: " + candidateMemorySize + " bytes");
            System.out.println("geneTreeCluster1 Memory: " + geneTreeMemorySize + " bytes");
            System.out.println("geneTreeCluster2 Memory: " + geneTreeMemorySize + " bytes");
            System.out.println("Total BitSet Memory: " + (2 * candidateMemorySize + 2 * geneTreeMemorySize) + " bytes");
            System.out.println("frequencies array size: " + (numGeneTreeBips * 4) + " bytes");
            System.out.println("weights array size: " + (numCandidates * 8) + " bytes");
            
            // Check for invalid parameters that could cause allocation errors
            if (numCandidates <= 0) {
                System.err.println("ERROR: numCandidates is " + numCandidates + " (must be > 0)");
                throw new RuntimeException("Invalid numCandidates: " + numCandidates);
            }
            if (numGeneTreeBips <= 0) {
                System.err.println("ERROR: numGeneTreeBips is " + numGeneTreeBips + " (must be > 0)");
                throw new RuntimeException("Invalid numGeneTreeBips: " + numGeneTreeBips);
            }
            if (bitsetSize <= 0) {
                System.err.println("ERROR: bitsetSize is " + bitsetSize + " (must be > 0)");
                throw new RuntimeException("Invalid bitsetSize: " + bitsetSize);
            }
            if (candidateMemorySize <= 0) {
                System.err.println("ERROR: candidateMemorySize is " + candidateMemorySize + " (must be > 0)");
                throw new RuntimeException("Invalid candidate memory size: " + candidateMemorySize);
            }
            if (geneTreeMemorySize <= 0) {
                System.err.println("ERROR: geneTreeMemorySize is " + geneTreeMemorySize + " (must be > 0)");
                throw new RuntimeException("Invalid gene tree memory size: " + geneTreeMemorySize);
            }
            
            System.out.println("==== PARAMETER VALIDATION PASSED ====");
            
            // Allocate arrays for CUDA with proper memory alignment
            System.out.println("==== ALLOCATING JNA STRUCTURES ====");
            System.out.println("Allocating cudaCandidates array of size: " + numCandidates);
            WeightCalcLib.Bipartition[] cudaCandidates = (WeightCalcLib.Bipartition[]) new WeightCalcLib.Bipartition().toArray(numCandidates);
            System.out.println("Successfully allocated cudaCandidates");
            
            System.out.println("Allocating cudaGeneTreeBips array of size: " + numGeneTreeBips);
            WeightCalcLib.Bipartition[] cudaGeneTreeBips = (WeightCalcLib.Bipartition[]) new WeightCalcLib.Bipartition().toArray(numGeneTreeBips);
            System.out.println("Successfully allocated cudaGeneTreeBips");
            
            System.out.println("Allocating frequencies array of size: " + numGeneTreeBips);
            int[] frequencies = new int[numGeneTreeBips];
            System.out.println("Successfully allocated frequencies array");
            
            System.out.println("Allocating weights array of size: " + numCandidates);
            double[] weights = new double[numCandidates];
            System.out.println("Successfully allocated weights array");
            
            // Allocate contiguous memory for all bit arrays
            System.out.println("==== ALLOCATING CONTIGUOUS MEMORY ====");
            System.out.println("Allocating candidateCluster1 Memory: " + candidateMemorySize + " bytes");
            Memory candidateCluster1 = new Memory(candidateMemorySize);
            System.out.println("Successfully allocated candidateCluster1");
            
            System.out.println("Allocating candidateCluster2 Memory: " + candidateMemorySize + " bytes");
            Memory candidateCluster2 = new Memory(candidateMemorySize);
            System.out.println("Successfully allocated candidateCluster2");
            
            System.out.println("Allocating geneTreeCluster1 Memory: " + geneTreeMemorySize + " bytes");
            Memory geneTreeCluster1 = new Memory(geneTreeMemorySize);
            System.out.println("Successfully allocated geneTreeCluster1");
            
            System.out.println("Allocating geneTreeCluster2 Memory: " + geneTreeMemorySize + " bytes");
            Memory geneTreeCluster2 = new Memory(geneTreeMemorySize);
            System.out.println("Successfully allocated geneTreeCluster2");
            
            System.out.println("==== ALL MEMORY ALLOCATIONS SUCCESSFUL ====");
            
            // Convert candidates
            System.out.println("==== CONVERTING CANDIDATES ====");
            for (int i = 0; i < numCandidates; i++) {
                STBipartition stb = candidates.get(i);
                if (stb == null) {
                    System.err.println("ERROR: Null candidate at index " + i);
                    throw new RuntimeException("Null candidate at index " + i);
                }
                if (stb.cluster1 == null) {
                    System.err.println("ERROR: Null cluster1 in candidate " + i);
                    throw new RuntimeException("Null cluster1 in candidate " + i);
                }
                if (stb.cluster2 == null) {
                    System.err.println("ERROR: Null cluster2 in candidate " + i);
                    throw new RuntimeException("Null cluster2 in candidate " + i);
                }
                
                WeightCalcLib.Bipartition bip = cudaCandidates[i];
                
                // Copy bit arrays to contiguous memory
                long[] cluster1Words = stb.cluster1.toLongArray();
                long[] cluster2Words = stb.cluster2.toLongArray();
                
                if (i < 3) { // Log first few for debugging
                    System.out.println("Candidate " + i + " cluster1 words: " + cluster1Words.length + 
                                     ", cluster2 words: " + cluster2Words.length);
                    System.out.println("Candidate " + i + " cluster1 cardinality: " + stb.cluster1.cardinality() + 
                                     ", cluster2 cardinality: " + stb.cluster2.cardinality());
                }
                
                candidateCluster1.write(i * bitsetSize * 8, cluster1Words, 0, cluster1Words.length);
                candidateCluster2.write(i * bitsetSize * 8, cluster2Words, 0, cluster2Words.length);
                
                // Set pointers to contiguous memory
                bip.cluster1 = candidateCluster1.share(i * bitsetSize * 8);
                bip.cluster2 = candidateCluster2.share(i * bitsetSize * 8);
                bip.bitsetSize = bitsetSize;
                bip.write();
            }
            
            System.out.println("==== CANDIDATES CONVERSION SUCCESSFUL ====");
            
            // Convert gene tree bipartitions
            System.out.println("==== CONVERTING GENE TREE BIPARTITIONS ====");
            int idx = 0;
            int totalFrequency = 0;
            for (Map.Entry<STBipartition, Integer> entry : geneTrees.stBipartitions.entrySet()) {
                STBipartition stb = entry.getKey();
                if (stb == null) {
                    System.err.println("ERROR: Null gene tree bipartition at index " + idx);
                    throw new RuntimeException("Null gene tree bipartition at index " + idx);
                }
                if (stb.cluster1 == null) {
                    System.err.println("ERROR: Null cluster1 in gene tree bipartition " + idx);
                    throw new RuntimeException("Null cluster1 in gene tree bipartition " + idx);
                }
                if (stb.cluster2 == null) {
                    System.err.println("ERROR: Null cluster2 in gene tree bipartition " + idx);
                    throw new RuntimeException("Null cluster2 in gene tree bipartition " + idx);
                }
                
                WeightCalcLib.Bipartition bip = cudaGeneTreeBips[idx];
                
                // Copy bit arrays to contiguous memory
                long[] cluster1Words = stb.cluster1.toLongArray();
                long[] cluster2Words = stb.cluster2.toLongArray();
                
                if (idx < 3) { // Log first few for debugging
                    System.out.println("GeneTree " + idx + " cluster1 words: " + cluster1Words.length + 
                                     ", cluster2 words: " + cluster2Words.length);
                    System.out.println("GeneTree " + idx + " cluster1 cardinality: " + stb.cluster1.cardinality() + 
                                     ", cluster2 cardinality: " + stb.cluster2.cardinality());
                    System.out.println("GeneTree " + idx + " frequency: " + entry.getValue());
                }
                
                geneTreeCluster1.write(idx * bitsetSize * 8, cluster1Words, 0, cluster1Words.length);
                geneTreeCluster2.write(idx * bitsetSize * 8, cluster2Words, 0, cluster2Words.length);
                
                // Set pointers to contiguous memory
                bip.cluster1 = geneTreeCluster1.share(idx * bitsetSize * 8);
                bip.cluster2 = geneTreeCluster2.share(idx * bitsetSize * 8);
                bip.bitsetSize = bitsetSize;
                bip.write();
                
                frequencies[idx] = entry.getValue();
                totalFrequency += entry.getValue();
                idx++;
            }
            
            System.out.println("==== GENE TREE CONVERSION SUCCESSFUL ====");
            System.out.println("Total frequency sum: " + totalFrequency);
            System.out.println("Processed " + idx + " unique gene tree bipartitions");
            
            // Launch CUDA kernel
            System.out.println("==== LAUNCHING CUDA KERNEL ====");
            System.out.println("Final kernel parameters:");
            System.out.println("  cudaCandidates: array of " + numCandidates + " elements");
            System.out.println("  cudaGeneTreeBips: array of " + numGeneTreeBips + " elements");
            System.out.println("  frequencies: array of " + numGeneTreeBips + " elements");
            System.out.println("  weights: array of " + numCandidates + " elements");
            System.out.println("  numCandidates: " + numCandidates);
            System.out.println("  numGeneTreeBips: " + numGeneTreeBips);
            System.out.println("  bitsetSize: " + bitsetSize);
            
            // Verify arrays are not null before kernel launch
            if (cudaCandidates == null) {
                System.err.println("ERROR: cudaCandidates is null!");
                throw new RuntimeException("cudaCandidates is null");
            }
            if (cudaGeneTreeBips == null) {
                System.err.println("ERROR: cudaGeneTreeBips is null!");
                throw new RuntimeException("cudaGeneTreeBips is null");
            }
            if (frequencies == null) {
                System.err.println("ERROR: frequencies is null!");
                throw new RuntimeException("frequencies is null");
            }
            if (weights == null) {
                System.err.println("ERROR: weights is null!");
                throw new RuntimeException("weights is null");
            }
            
            System.out.println("About to call CUDA kernel...");
            
            WeightCalcLib.INSTANCE.launchWeightCalculation(
                cudaCandidates,
                cudaGeneTreeBips,
                frequencies,
                weights,
                numCandidates,
                numGeneTreeBips,
                bitsetSize
            );
            
            System.out.println("==== CUDA KERNEL COMPLETED SUCCESSFULLY ====");
            
            // Convert results back to Java Map
            System.out.println("==== CONVERTING RESULTS BACK TO JAVA ====");
            Map<STBipartition, Double> result = new HashMap<>();
            double minWeight = Double.MAX_VALUE;
            double maxWeight = Double.MIN_VALUE;
            double totalWeight = 0.0;
            int zeroWeights = 0;
            int negativeWeights = 0;
            
            for (int i = 0; i < numCandidates; i++) {
                double weight = weights[i];
                result.put(candidates.get(i), weight);
                
                // Collect statistics
                minWeight = Math.min(minWeight, weight);
                maxWeight = Math.max(maxWeight, weight);
                totalWeight += weight;
                if (weight == 0.0) zeroWeights++;
                if (weight < 0.0) negativeWeights++;
            }
            
            // Print comprehensive result statistics
            System.out.println("==== WEIGHT CALCULATION STATISTICS ====");
            System.out.println("Total weights computed: " + numCandidates);
            System.out.println("Min weight: " + minWeight);
            System.out.println("Max weight: " + maxWeight);
            System.out.println("Average weight: " + (totalWeight / numCandidates));
            System.out.println("Zero weights: " + zeroWeights);
            System.out.println("Negative weights: " + negativeWeights);
            
            // Print some sample weights for debugging
            System.out.println("==== SAMPLE WEIGHTS ====");
            int count = 0;
            for (Map.Entry<STBipartition, Double> entry : result.entrySet()) {
                if (count < 5) {
                    System.out.println("Weight " + count + ": " + entry.getValue());
                    count++;
                } else {
                    break;
                }
            }
            
            System.out.println("==== GPU WEIGHT CALCULATION COMPLETED SUCCESSFULLY ====");
            return result;
            
        } catch (Exception e) {
            System.err.println("==== GPU COMPUTATION FAILED ====");
            System.err.println("Error message: " + e.getMessage());
            System.err.println("Exception type: " + e.getClass().getSimpleName());
            System.err.println("Stack trace:");
            e.printStackTrace();
            System.err.println("==== FALLING BACK TO CPU COMPUTATION ====");
            return calculateWeightsMultiThread(candidates);
        }
    }
    
    private Pointer convertBitSetToPointer(BitSet bitSet) {
        // Convert Java BitSet to native memory
        long[] words = bitSet.toLongArray();
        Pointer ptr = new Memory(words.length * 8); // 8 bytes per long
        ptr.write(0, words, 0, words.length);
        return ptr;
    }
} 