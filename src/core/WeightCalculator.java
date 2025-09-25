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
            // Use long for all size calculations to prevent overflow
            long numCandidatesLong = candidates.size();
            long numGeneTreeBipsLong = geneTrees.stBipartitions.size();
            // Calculate bitset size safely to prevent overflow
            long realTaxaCountLong = geneTrees.realTaxaCount;
            long bitsetSizeLong = (realTaxaCountLong + 63L) / 64L; // Size in 64-bit words
            
            // Extra safety: ensure realTaxaCount itself is reasonable
            if (realTaxaCountLong < 0) {
                throw new RuntimeException("Invalid realTaxaCount: " + realTaxaCountLong + " (cannot be negative)");
            }
            if (realTaxaCountLong > 1000000000L) { // 1 billion taxa limit
                throw new RuntimeException("Too many taxa: " + realTaxaCountLong + " exceeds reasonable limit of 1 billion");
            }
            
            // Validate that sizes fit in int for array indexing (Java arrays are int-indexed)
            if (numCandidatesLong > Integer.MAX_VALUE) {
                throw new RuntimeException("Too many candidates: " + numCandidatesLong + " exceeds Integer.MAX_VALUE");
            }
            if (numGeneTreeBipsLong > Integer.MAX_VALUE) {
                throw new RuntimeException("Too many gene tree bipartitions: " + numGeneTreeBipsLong + " exceeds Integer.MAX_VALUE");
            }
            if (bitsetSizeLong > Integer.MAX_VALUE) {
                throw new RuntimeException("Bitset size too large: " + bitsetSizeLong + " exceeds Integer.MAX_VALUE");
            }
            
            // Safe to cast to int for array indexing after validation
            int numCandidates = (int) numCandidatesLong;
            int numGeneTreeBips = (int) numGeneTreeBipsLong;
            int bitsetSize = (int) bitsetSizeLong;
            
            System.out.println("==== CUDA ALLOCATION PARAMETERS ====");
            System.out.println("numCandidates: " + numCandidates);
            System.out.println("numGeneTreeBips: " + numGeneTreeBips);
            System.out.println("bitsetSize (in 64-bit words): " + bitsetSize);
            System.out.println("geneTrees.realTaxaCount: " + geneTrees.realTaxaCount);
            
            // Calculate memory allocation sizes using long arithmetic to prevent overflow
            // Use the original long values for all calculations
            long candidateMemorySize = numCandidatesLong * bitsetSizeLong * 8L;
            long geneTreeMemorySize = numGeneTreeBipsLong * bitsetSizeLong * 8L;
            long frequenciesSize = numGeneTreeBipsLong * 4L;
            long weightsSize = numCandidatesLong * 8L;
            long totalMemory = 2L * candidateMemorySize + 2L * geneTreeMemorySize + frequenciesSize + weightsSize;
            
            System.out.println("==== MEMORY ALLOCATION SIZES (using long arithmetic) ====");
            System.out.println("candidateCluster1 Memory: " + candidateMemorySize + " bytes (" + (candidateMemorySize / (1024*1024)) + " MB)");
            System.out.println("candidateCluster2 Memory: " + candidateMemorySize + " bytes (" + (candidateMemorySize / (1024*1024)) + " MB)");
            System.out.println("geneTreeCluster1 Memory: " + geneTreeMemorySize + " bytes (" + (geneTreeMemorySize / (1024*1024)) + " MB)");
            System.out.println("geneTreeCluster2 Memory: " + geneTreeMemorySize + " bytes (" + (geneTreeMemorySize / (1024*1024)) + " MB)");
            System.out.println("frequencies array size: " + frequenciesSize + " bytes");
            System.out.println("weights array size: " + weightsSize + " bytes");
            System.out.println("Total Memory Required: " + totalMemory + " bytes (" + (totalMemory / (1024*1024)) + " MB)");
            
            // Check for potential integer overflow in calculations (demonstrate what would happen with int arithmetic)
            System.out.println("==== OVERFLOW DETECTION ====");
            long candidateCalcInt = (long) numCandidates * bitsetSize * 8; // Still potential int overflow in multiplication
            long geneTreeCalcInt = (long) numGeneTreeBips * bitsetSize * 8;  // Still potential int overflow in multiplication
            System.out.println("candidateMemorySize (mixed int arithmetic): " + candidateCalcInt + " (overflow if negative!)");
            System.out.println("geneTreeMemorySize (mixed int arithmetic): " + geneTreeCalcInt + " (overflow if negative!)");
            System.out.println("candidateMemorySize (pure long arithmetic): " + candidateMemorySize + " (always correct)");
            System.out.println("geneTreeMemorySize (pure long arithmetic): " + geneTreeMemorySize + " (always correct)");
            
            if (candidateCalcInt != candidateMemorySize) {
                System.err.println("WARNING: Integer overflow would occur in candidate memory calculation!");
                System.err.println("  Mixed arithmetic result: " + candidateCalcInt + ", Pure long result: " + candidateMemorySize);
            }
            if (geneTreeCalcInt != geneTreeMemorySize) {
                System.err.println("WARNING: Integer overflow would occur in gene tree memory calculation!");
                System.err.println("  Mixed arithmetic result: " + geneTreeCalcInt + ", Pure long result: " + geneTreeMemorySize);
            }
            
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
            
            // Check for memory sizes that exceed JNA Memory constructor limits
            // JNA Memory constructor takes long but may have implementation limits
            long maxReasonableSize = 8L * 1024 * 1024 * 1024; // 8 GB limit
            if (candidateMemorySize > maxReasonableSize) {
                System.err.println("ERROR: candidateMemorySize (" + candidateMemorySize + 
                                 " bytes = " + (candidateMemorySize / (1024*1024)) + " MB) exceeds reasonable limit");
                throw new RuntimeException("Candidate memory size too large: " + candidateMemorySize + " bytes");
            }
            if (geneTreeMemorySize > maxReasonableSize) {
                System.err.println("ERROR: geneTreeMemorySize (" + geneTreeMemorySize + 
                                 " bytes = " + (geneTreeMemorySize / (1024*1024)) + " MB) exceeds reasonable limit");
                throw new RuntimeException("Gene tree memory size too large: " + geneTreeMemorySize + " bytes");
            }
            if (totalMemory > 16L * 1024 * 1024 * 1024) { // 16 GB total limit
                System.err.println("ERROR: Total memory requirement (" + totalMemory + 
                                 " bytes = " + (totalMemory / (1024*1024)) + " MB) exceeds reasonable limit");
                throw new RuntimeException("Total memory requirement too large: " + totalMemory + " bytes");
            }
            
            // Calculate when offset overflow would occur with int arithmetic
            System.out.println("==== OVERFLOW THRESHOLD ANALYSIS ====");
            long maxSafeIndex = Integer.MAX_VALUE / (bitsetSizeLong * 8L);
            System.out.println("Max safe candidate/gene tree index before offset overflow: " + maxSafeIndex);
            System.out.println("Actual numCandidates: " + numCandidatesLong + (numCandidatesLong > maxSafeIndex ? " [WOULD OVERFLOW WITH INT ARITHMETIC!]" : " [safe]"));
            System.out.println("Actual numGeneTreeBips: " + numGeneTreeBipsLong + (numGeneTreeBipsLong > maxSafeIndex ? " [WOULD OVERFLOW WITH INT ARITHMETIC!]" : " [safe]"));
            
            if (numCandidatesLong > maxSafeIndex) {
                System.err.println("INFO: numCandidates (" + numCandidatesLong + ") exceeds int-safe limit (" + maxSafeIndex + ")");
                System.err.println("      But we're using long arithmetic, so no overflow will occur!");
            }
            if (numGeneTreeBipsLong > maxSafeIndex) {
                System.err.println("INFO: numGeneTreeBips (" + numGeneTreeBipsLong + ") exceeds int-safe limit (" + maxSafeIndex + ")");
                System.err.println("      But we're using long arithmetic, so no overflow will occur!");
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
                
                // Use long arithmetic to prevent overflow in memory offset calculations
                long candidateOffset = (long) i * bitsetSizeLong * 8L;
                
                // Debug overflow detection for offset calculation
                if (i < 3 || candidateOffset != (long)(i * bitsetSize * 8)) {
                    int offsetInt = i * bitsetSize * 8; // int arithmetic - may overflow
                    System.out.println("Candidate " + i + " offset: long=" + candidateOffset + 
                                     ", int=" + offsetInt + (offsetInt != candidateOffset ? " [OVERFLOW!]" : " [OK]"));
                    
                    if (offsetInt != candidateOffset) {
                        System.err.println("INFO: Offset overflow would occur at candidate " + i + 
                                         " with int arithmetic (but we're using long, so it's safe)");
                        System.err.println("  Calculation: " + i + " * " + bitsetSizeLong + " * 8 = " + 
                                         candidateOffset + " (long, correct) vs " + offsetInt + " (int, overflowed)");
                    }
                }
                
                candidateCluster1.write(candidateOffset, cluster1Words, 0, cluster1Words.length);
                candidateCluster2.write(candidateOffset, cluster2Words, 0, cluster2Words.length);
                
                // Set pointers to contiguous memory
                bip.cluster1 = candidateCluster1.share(candidateOffset);
                bip.cluster2 = candidateCluster2.share(candidateOffset);
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
                
                // Use long arithmetic to prevent overflow in memory offset calculations
                long geneTreeOffset = (long) idx * bitsetSizeLong * 8L;
                
                // Debug overflow detection for offset calculation
                if (idx < 3 || geneTreeOffset != (long)(idx * bitsetSize * 8)) {
                    int offsetInt = idx * bitsetSize * 8; // int arithmetic - may overflow
                    System.out.println("GeneTree " + idx + " offset: long=" + geneTreeOffset + 
                                     ", int=" + offsetInt + (offsetInt != geneTreeOffset ? " [OVERFLOW!]" : " [OK]"));
                    
                    if (offsetInt != geneTreeOffset) {
                        System.err.println("INFO: Offset overflow would occur at gene tree " + idx + 
                                         " with int arithmetic (but we're using long, so it's safe)");
                        System.err.println("  Calculation: " + idx + " * " + bitsetSizeLong + " * 8 = " + 
                                         geneTreeOffset + " (long, correct) vs " + offsetInt + " (int, overflowed)");
                    }
                }
                
                geneTreeCluster1.write(geneTreeOffset, cluster1Words, 0, cluster1Words.length);
                geneTreeCluster2.write(geneTreeOffset, cluster2Words, 0, cluster2Words.length);
                
                // Set pointers to contiguous memory
                bip.cluster1 = geneTreeCluster1.share(geneTreeOffset);
                bip.cluster2 = geneTreeCluster2.share(geneTreeOffset);
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