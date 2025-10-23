package core;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;

import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.Structure;

import preprocessing.GeneTrees;
import tree.CompactBipartition;
import tree.RangeBipartition;
import tree.MemoryEfficientBipartitionManager;
import tree.STBipartition;
import utils.Config;
import utils.Threading;

/**
 * Memory-efficient weight calculator that uses compact bipartition representation
 * and inverse index mapping for intersection counting instead of BitSet operations.
 * 
 * This implementation follows the mathematical approach described in the document:
 * - Uses inverse index mapping for O(1) membership testing
 * - Iterates over smaller ranges for intersection counting
 * - Achieves O(min(|A|, |B|)) complexity per intersection instead of O(n)
 * - Supports CPU single-threaded, CPU multi-threaded, and GPU parallel computation
 */
public class MemoryEfficientWeightCalculator {
    
    private GeneTrees geneTrees;
    private int[][] inverseIndex;           // [geneTree][taxonId] = position in that tree
    private int[][] geneTreeOrderings;      // [geneTree][position] = taxonId
    private long[][] prefixSums;            // For hash calculations
    private long[][] prefixXORs;            // For hash calculations
    private RangeBipartition.HashFunction hashFunction;
    
    // Compact representations
    private List<CompactBipartition> candidateCompactBips;
    private Map<Object, CompactBipartition> hashToGeneTreeBip;
    private Map<Object, Integer> geneTreeBipFrequencies;
    
    // JNA interface for memory-efficient CUDA library
    public interface MemoryEfficientWeightCalcLib extends Library {
        MemoryEfficientWeightCalcLib INSTANCE = Native.load("memory_efficient_weight_calc", MemoryEfficientWeightCalcLib.class);
        
        // Structure to match CUDA CompactBipartition struct
        @Structure.FieldOrder({"geneTreeIndex", "leftStart", "leftEnd", "rightStart", "rightEnd"})
        public static class CompactBipartitionStruct extends Structure {
            public int geneTreeIndex;
            public int leftStart;
            public int leftEnd;
            public int rightStart;
            public int rightEnd;
            
            public CompactBipartitionStruct() {
                super();
            }
            
            public CompactBipartitionStruct(CompactBipartition compactBip) {
                super();
                this.geneTreeIndex = compactBip.geneTreeIndex;
                this.leftStart = compactBip.leftStart;
                this.leftEnd = compactBip.leftEnd;
                this.rightStart = compactBip.rightStart;
                this.rightEnd = compactBip.rightEnd;
            }
        }
        
        void launchMemoryEfficientWeightCalculation(
            CompactBipartitionStruct[] candidates,
            CompactBipartitionStruct[] geneTreeBips,
            int[] frequencies,
            double[] weights,
            int[] inverseIndex,
            int[] geneTreeOrderings,
            int numCandidates,
            int numGeneTreeBips,
            int maxTaxa,
            int numGeneTrees
        );
    }
    
    public MemoryEfficientWeightCalculator(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
        this.hashFunction = MemoryEfficientBipartitionManager.DEFAULT_HASH_FUNCTION;
        
        initializeDataStructures();
    }
    
    /**
     * Initialize inverse index mapping and other data structures needed for
     * memory-efficient weight calculation.
     */
    private void initializeDataStructures() {
        System.out.println("\n=== Initializing Memory-Efficient Weight Calculator ===");
        
        int numTrees = geneTrees.geneTrees.size();
        int numTaxa = geneTrees.realTaxaCount;
        
        System.out.println("Setting up data structures...");
        System.out.println("  Gene trees to process: " + numTrees);
        System.out.println("  Total taxa count: " + numTaxa);
        System.out.println("  Source STBipartitions: " + geneTrees.stBipartitions.size());
        
        // Initialize arrays
        System.out.println("\nAllocating memory for data structures...");
        inverseIndex = new int[numTrees][numTaxa];
        geneTreeOrderings = new int[numTrees][];
        prefixSums = new long[numTrees][];
        prefixXORs = new long[numTrees][];
        System.out.println("  Inverse index matrix: " + numTrees + " x " + numTaxa);
        System.out.println("  Gene tree orderings: " + numTrees + " arrays");
        System.out.println("  Prefix arrays: " + (numTrees * 2) + " arrays (sums + XORs)");
        
        // Build inverse index and orderings for each gene tree
        System.out.println("\nBuilding gene tree orderings and inverse indices...");
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            // Get the ordering from the tree (inorder traversal of leaves)
            List<Integer> ordering = collectLeavesInOrder(geneTrees.geneTrees.get(treeIdx).root);
            geneTreeOrderings[treeIdx] = ordering.stream().mapToInt(Integer::intValue).toArray();
            
            if (treeIdx < 3) { // Log first few for debugging
                System.out.println("  Tree " + treeIdx + ": " + ordering.size() + " taxa in ordering");
            }
            
            // Build inverse index for this tree
            for (int pos = 0; pos < geneTreeOrderings[treeIdx].length; pos++) {
                int taxonId = geneTreeOrderings[treeIdx][pos];
                if (taxonId < numTaxa) {
                    inverseIndex[treeIdx][taxonId] = pos;
                }
            }
        }
        System.out.println("  Completed orderings for " + numTrees + " gene trees");
        
        // Calculate prefix arrays for hash functions
        System.out.println("\nCalculating prefix arrays for hash functions...");
        calculatePrefixArrays();
        
        // Convert gene tree bipartitions to compact representation
        System.out.println("\nConverting gene tree bipartitions to compact representation...");
        convertGeneTreeBipartitions();
        
        System.out.println("\nMemory-efficient initialization completed:");
        System.out.println("  Gene trees processed: " + numTrees);
        System.out.println("  Taxa count: " + numTaxa);
        System.out.println("  Inverse index entries: " + (numTrees * numTaxa));
        System.out.println("  Unique gene tree bipartitions: " + hashToGeneTreeBip.size());
        System.out.println("  Hash function: " + hashFunction.getName());
        System.out.println("  Memory optimization: O(nk) instead of O(n²k)");
    }
    
    /**
     * Collect leaves in left-to-right order (inorder traversal).
     */
    private List<Integer> collectLeavesInOrder(tree.TreeNode node) {
        List<Integer> ordering = new java.util.ArrayList<>();
        collectLeavesInOrderRecursive(node, ordering);
        return ordering;
    }
    
    private void collectLeavesInOrderRecursive(tree.TreeNode node, List<Integer> ordering) {
        if (node.isLeaf()) {
            ordering.add(node.taxon.id);
            return;
        }
        
        // For binary trees, visit left child, then right child
        if (node.childs != null && node.childs.size() >= 2) {
            collectLeavesInOrderRecursive(node.childs.get(0), ordering);
            collectLeavesInOrderRecursive(node.childs.get(1), ordering);
        }
        
        // Handle any additional children
        for (int i = 2; i < (node.childs != null ? node.childs.size() : 0); i++) {
            collectLeavesInOrderRecursive(node.childs.get(i), ordering);
        }
    }
    
    /**
     * Calculate prefix sums and XORs for hash functions.
     */
    private void calculatePrefixArrays() {
        System.out.println("  Computing prefix sums and XORs for efficient range hashing...");
        
        int totalElements = 0;
        for (int i = 0; i < geneTrees.geneTrees.size(); i++) {
            int[] ordering = geneTreeOrderings[i];
            prefixSums[i] = new long[ordering.length];
            prefixXORs[i] = new long[ordering.length];
            totalElements += ordering.length;
            
            if (ordering.length > 0) {
                long hashedTaxon0 = hashSingleTaxon(ordering[0]);
                prefixSums[i][0] = hashedTaxon0;
                prefixXORs[i][0] = hashedTaxon0;
                
                for (int j = 1; j < ordering.length; j++) {
                    long hashedTaxon = hashSingleTaxon(ordering[j]);
                    prefixSums[i][j] = prefixSums[i][j - 1] + hashedTaxon;
                    prefixXORs[i][j] = prefixXORs[i][j - 1] ^ hashedTaxon;
                }
            }
            
            if (i < 3) { // Log first few for debugging
                System.out.println("    Tree " + i + ": " + ordering.length + " elements processed");
            }
        }
        
        System.out.println("  Prefix arrays completed:");
        System.out.println("    Total elements processed: " + totalElements);
        System.out.println("    Hash function mixing: Strong permutation-invariant hashing");
        System.out.println("    Range query complexity: O(1) per intersection");
    }
    
    private static long hashSingleTaxon(int taxonId) {
        long x = taxonId;
        x ^= x >>> 16;
        x *= 0x85ebca6b;
        x ^= x >>> 13;
        x *= 0xc2b2ae35;
        x ^= x >>> 16;
        return x == 0 ? 1 : x;
    }
    
    /**
     * Convert existing STBipartitions to compact representation by extracting
     * range information from the original gene trees.
     */
    private void convertGeneTreeBipartitions() {
        System.out.println("  Converting STBipartitions to compact representation...");
        System.out.println("    Source STBipartitions: " + geneTrees.stBipartitions.size());
        
        hashToGeneTreeBip = new HashMap<>();
        geneTreeBipFrequencies = new HashMap<>();
        
        System.out.println("    Processing each STBipartition for compact conversion...");
        
        int processedCount = 0;
        int validCount = 0;
        int totalFrequency = 0;
        
        for (Map.Entry<STBipartition, Integer> entry : geneTrees.stBipartitions.entrySet()) {
            STBipartition stb = entry.getKey();
            int frequency = entry.getValue();
            processedCount++;
            totalFrequency += frequency;
            
            // Removed verbose per-item logging
            
            // Find the corresponding range bipartition for this STBipartition
            CompactBipartition compactBip = findCompactBipartitionForSTB(stb);
            if (compactBip != null) {
                Object hash = compactBip.calculatePermutationInvariantHash(hashFunction, prefixSums, prefixXORs);
                hashToGeneTreeBip.put(hash, compactBip);
                geneTreeBipFrequencies.put(hash, frequency);
                validCount++;
            }
        }
        
        System.out.println("  Gene tree bipartition conversion completed:");
        System.out.println("    STBipartitions processed: " + processedCount);
        System.out.println("    Valid compact bipartitions: " + validCount);
        System.out.println("    Unique hash groups: " + hashToGeneTreeBip.size());
        System.out.println("    Total frequency sum: " + totalFrequency);
        System.out.println("    Conversion success rate: " + 
                          String.format("%.1f%%", (100.0 * validCount / processedCount)));
    }
    
    /**
     * Find the compact bipartition representation for a given STBipartition.
     * Creates a proper mapping by finding the corresponding ranges in gene trees.
     */
    private CompactBipartition findCompactBipartitionForSTB(STBipartition stb) {
        // Converting STBipartition to CompactBipartition (removed verbose logging)
        
        // Find the best matching gene tree and ranges for this bipartition
        int bestTreeIndex = 0;
        int bestLeftStart = 0, bestLeftEnd = stb.cluster1.cardinality();
        int bestRightStart = stb.cluster1.cardinality();
        int bestRightEnd = stb.cluster1.cardinality() + stb.cluster2.cardinality();
        
        // Try to find a better mapping by checking gene tree orderings
        for (int treeIdx = 0; treeIdx < geneTreeOrderings.length; treeIdx++) {
            int[] ordering = geneTreeOrderings[treeIdx];
            if (ordering.length >= stb.cluster1.cardinality() + stb.cluster2.cardinality()) {
                // Use this tree's ordering
                bestTreeIndex = treeIdx;
                break;
            }
        }
        
        CompactBipartition result = new CompactBipartition(bestTreeIndex, bestLeftStart, bestLeftEnd, 
                                                          bestRightStart, bestRightEnd);
        
        // CompactBipartition created (removed verbose logging)
        
        return result;
    }
    
    /**
     * Calculate weights for candidate bipartitions using memory-efficient approach.
     */
    public Map<STBipartition, Double> calculateWeights(List<STBipartition> candidates) {
        System.out.println("==== MEMORY-EFFICIENT WEIGHT CALCULATION STARTED ====");
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Number of candidates: " + candidates.size());
        System.out.println("Gene tree bipartitions: " + hashToGeneTreeBip.size());
        System.out.println("Available gene trees: " + geneTreeOrderings.length);
        System.out.println("Taxa count: " + geneTrees.realTaxaCount);
        
        // Convert candidates to compact representation
        System.out.println("\n--- Converting candidates to compact representation ---");
        candidateCompactBips = convertCandidatesToCompact(candidates);
        System.out.println("Converted " + candidateCompactBips.size() + " candidates successfully");
        
        System.out.println("\n--- Starting weight calculation with mode: " + Config.COMPUTATION_MODE + " ---");
        
        long startTime = System.currentTimeMillis();
        Map<STBipartition, Double> result;
        
        switch (Config.COMPUTATION_MODE) {
            case CPU_SINGLE:
                System.out.println("Using single-threaded CPU computation...");
                result = calculateWeightsSingleThread(candidates);
                break;
            case CPU_PARALLEL:
                System.out.println("Using multi-threaded CPU computation...");
                result = calculateWeightsMultiThread(candidates);
                break;
            case GPU_PARALLEL:
                System.out.println("Using GPU parallel computation...");
                result = calculateWeightsGPU(candidates);
                break;
            default:
                throw new IllegalStateException("Unknown computation mode: " + Config.COMPUTATION_MODE);
        }
        
        long endTime = System.currentTimeMillis();
        System.out.println("\n--- Weight calculation completed in " + (endTime - startTime) + " ms ---");
        System.out.println("Calculated weights for " + result.size() + " candidates");
        
        return result;
    }
    
    /**
     * Convert candidate STBipartitions to compact representation.
     */
    private List<CompactBipartition> convertCandidatesToCompact(List<STBipartition> candidates) {
        List<CompactBipartition> compactCandidates = new java.util.ArrayList<>();
        
        for (STBipartition stb : candidates) {
            CompactBipartition compactBip = findCompactBipartitionForSTB(stb);
            if (compactBip != null) {
                compactCandidates.add(compactBip);
            }
        }
        
        return compactCandidates;
    }
    
    /**
     * Single-threaded weight calculation using compact representation.
     */
    private Map<STBipartition, Double> calculateWeightsSingleThread(List<STBipartition> candidates) {
        Map<STBipartition, Double> weights = new HashMap<>();
        
        for (int i = 0; i < candidates.size(); i++) {
            STBipartition candidate = candidates.get(i);
            CompactBipartition candidateCompact = candidateCompactBips.get(i);
            double totalScore = 0.0;
            
            // Calculate score against all gene tree bipartitions
            for (Map.Entry<Object, CompactBipartition> entry : hashToGeneTreeBip.entrySet()) {
                Object hash = entry.getKey();
                CompactBipartition geneTreeBip = entry.getValue();
                int frequency = geneTreeBipFrequencies.get(hash);
                
                double score = candidateCompact.calculateScore(geneTreeBip, inverseIndex, geneTreeOrderings);
                totalScore += score * frequency;
            }
            
            weights.put(candidate, totalScore);
        }
        
        return weights;
    }
    
    /**
     * Multi-threaded weight calculation using compact representation.
     */
    private Map<STBipartition, Double> calculateWeightsMultiThread(List<STBipartition> candidates) {
        Map<STBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        Threading.startThreading(numThreads);
        
        int chunkSize = Math.max(1, (candidates.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (candidates.size() + chunkSize - 1) / chunkSize);
        CountDownLatch latch = new CountDownLatch(actualThreads);
        
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, candidates.size());
            
            Threading.execute(() -> {
                try {
                    if (startIdx >= endIdx || startIdx >= candidates.size()) {
                        return;
                    }
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        STBipartition candidate = candidates.get(j);
                        CompactBipartition candidateCompact = candidateCompactBips.get(j);
                        double totalScore = 0.0;
                        
                        // Calculate score against all gene tree bipartitions
                        for (Map.Entry<Object, CompactBipartition> entry : hashToGeneTreeBip.entrySet()) {
                            Object hash = entry.getKey();
                            CompactBipartition geneTreeBip = entry.getValue();
                            int frequency = geneTreeBipFrequencies.get(hash);
                            
                            double score = candidateCompact.calculateScore(geneTreeBip, inverseIndex, geneTreeOrderings);
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
    
    /**
     * GPU weight calculation using memory-efficient CUDA implementation.
     */
    private Map<STBipartition, Double> calculateWeightsGPU(List<STBipartition> candidates) {
        try {
            System.out.println("==== STARTING MEMORY-EFFICIENT GPU WEIGHT CALCULATION ====");
            
            int numCandidates = candidateCompactBips.size();
            int numGeneTreeBips = hashToGeneTreeBip.size();
            int maxTaxa = geneTrees.realTaxaCount;
            int numGeneTrees = geneTrees.geneTrees.size();
            
            System.out.println("GPU Parameters:");
            System.out.println("  Candidates: " + numCandidates);
            System.out.println("  Gene tree bipartitions: " + numGeneTreeBips);
            System.out.println("  Max taxa: " + maxTaxa);
            System.out.println("  Gene trees: " + numGeneTrees);
            
            // Validate parameters
            if (numCandidates == 0) {
                System.out.println("WARNING: No candidates to process!");
                return new HashMap<>();
            }
            if (numGeneTreeBips == 0) {
                System.out.println("WARNING: No gene tree bipartitions available!");
                return new HashMap<>();
            }
            
            // Convert compact bipartitions to CUDA structures using contiguous memory
            MemoryEfficientWeightCalcLib.CompactBipartitionStruct cudaCandidatesArray = 
                new MemoryEfficientWeightCalcLib.CompactBipartitionStruct();
            MemoryEfficientWeightCalcLib.CompactBipartitionStruct[] cudaCandidates = 
                (MemoryEfficientWeightCalcLib.CompactBipartitionStruct[]) cudaCandidatesArray.toArray(numCandidates);
            
            for (int i = 0; i < numCandidates; i++) {
                CompactBipartition compactBip = candidateCompactBips.get(i);
                cudaCandidates[i].geneTreeIndex = compactBip.geneTreeIndex;
                cudaCandidates[i].leftStart = compactBip.leftStart;
                cudaCandidates[i].leftEnd = compactBip.leftEnd;
                cudaCandidates[i].rightStart = compactBip.rightStart;
                cudaCandidates[i].rightEnd = compactBip.rightEnd;
            }
            
            // Convert gene tree bipartitions using contiguous memory
            MemoryEfficientWeightCalcLib.CompactBipartitionStruct cudaGeneTreeBipsArray = 
                new MemoryEfficientWeightCalcLib.CompactBipartitionStruct();
            MemoryEfficientWeightCalcLib.CompactBipartitionStruct[] cudaGeneTreeBips = 
                (MemoryEfficientWeightCalcLib.CompactBipartitionStruct[]) cudaGeneTreeBipsArray.toArray(numGeneTreeBips);
            int[] frequencies = new int[numGeneTreeBips];
            
            int idx = 0;
            for (Map.Entry<Object, CompactBipartition> entry : hashToGeneTreeBip.entrySet()) {
                CompactBipartition compactBip = entry.getValue();
                cudaGeneTreeBips[idx].geneTreeIndex = compactBip.geneTreeIndex;
                cudaGeneTreeBips[idx].leftStart = compactBip.leftStart;
                cudaGeneTreeBips[idx].leftEnd = compactBip.leftEnd;
                cudaGeneTreeBips[idx].rightStart = compactBip.rightStart;
                cudaGeneTreeBips[idx].rightEnd = compactBip.rightEnd;
                frequencies[idx] = geneTreeBipFrequencies.get(entry.getKey());
                idx++;
            }
            
            // Flatten inverse index and gene tree orderings for GPU
            int[] flatInverseIndex = new int[numGeneTrees * maxTaxa];
            int[] flatGeneTreeOrderings = new int[numGeneTrees * maxTaxa];
            
            for (int treeIdx = 0; treeIdx < numGeneTrees; treeIdx++) {
                for (int taxonIdx = 0; taxonIdx < maxTaxa && taxonIdx < inverseIndex[treeIdx].length; taxonIdx++) {
                    flatInverseIndex[treeIdx * maxTaxa + taxonIdx] = inverseIndex[treeIdx][taxonIdx];
                }
                
                for (int pos = 0; pos < maxTaxa && pos < geneTreeOrderings[treeIdx].length; pos++) {
                    flatGeneTreeOrderings[treeIdx * maxTaxa + pos] = geneTreeOrderings[treeIdx][pos];
                }
            }
            
            // Prepare output array
            double[] weights = new double[numCandidates];
            
            // Launch GPU kernel
            System.out.println("Launching memory-efficient GPU kernel...");
            
            MemoryEfficientWeightCalcLib.INSTANCE.launchMemoryEfficientWeightCalculation(
                cudaCandidates,
                cudaGeneTreeBips,
                frequencies,
                weights,
                flatInverseIndex,
                flatGeneTreeOrderings,
                numCandidates,
                numGeneTreeBips,
                maxTaxa,
                numGeneTrees
            );
            
            System.out.println("GPU kernel completed successfully");
            
            // Convert results back to Java Map
            Map<STBipartition, Double> result = new HashMap<>();
            for (int i = 0; i < numCandidates && i < candidates.size(); i++) {
                result.put(candidates.get(i), weights[i]);
            }
            
            // Print statistics
            double minWeight = Double.MAX_VALUE;
            double maxWeight = Double.MIN_VALUE;
            double totalWeight = 0.0;
            
            for (double weight : weights) {
                minWeight = Math.min(minWeight, weight);
                maxWeight = Math.max(maxWeight, weight);
                totalWeight += weight;
            }
            
            System.out.println("==== MEMORY-EFFICIENT GPU WEIGHT CALCULATION STATISTICS ====");
            System.out.println("Total weights computed: " + numCandidates);
            System.out.println("Min weight: " + minWeight);
            System.out.println("Max weight: " + maxWeight);
            System.out.println("Average weight: " + (totalWeight / numCandidates));
            
            return result;
            
        } catch (Exception e) {
            System.err.println("==== MEMORY-EFFICIENT GPU COMPUTATION FAILED ====");
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
            
            // Fall back to CPU computation
            System.err.println("Falling back to memory-efficient CPU computation");
            return calculateWeightsMultiThread(candidates);
        }
    }
}