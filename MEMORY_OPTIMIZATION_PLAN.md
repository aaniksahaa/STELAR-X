# Memory-Optimized Weight Calculation and DP Implementation Plan

## Overview

This document outlines the implementation plan for optimizing the STELAR-MP phylogeny reconstruction algorithm to avoid memory-intensive BitSet expansions in both the weight calculation and dynamic programming (DP) phases. The optimization builds upon the existing memory-efficient bipartition representation using integer tuples and extends it to the remaining computational phases.

## Current System Analysis

### Current Implementation Status
1. **✅ Memory-Efficient Bipartition Generation**: Already implemented using `RangeBipartition` (integer tuples) and `MemoryEfficientBipartitionManager` with double hashing
2. **❌ Weight Calculation**: Currently uses BitSet expansion - converts all bipartitions to full BitSets for intersection calculations (memory-intensive)
3. **❌ DP Algorithm**: Uses BitSet clusters as HashMap keys, causing memory overhead
4. **❌ GPU Implementation**: Transfers full BitSet arrays to CUDA, requiring O(n²k) memory

### Problems to Solve
1. Weight calculation still expands to BitSets for intersection computation
2. DP uses BitSet clusters as keys instead of compact hash representations  
3. GPU weight calculation transfers massive BitSet arrays
4. Memory usage scales as O(n²k) instead of the theoretical O(nk)

## Proposed Solution Architecture

### Phase 1: Memory-Optimized Weight Calculation

#### 1.1 Create New Compact Data Structures

**`Cluster` Class** (`src/tree/Cluster.java`)
- Represents a cluster as `(geneTreeIndex, startPos, endPos)` with hash support
- Uses double hashing for collision resistance
- Isolated from existing `RangeBipartition` class but uses same hashing scheme

**`HashUtils` Class** (`src/utils/HashUtils.java`)  
- Shared hashing utilities using the same double-hashing scheme as `RangeBipartition`
- Provides cluster-specific hash computation
- Returns `HashPair` (sum hash, xor hash) for double hashing

**`InverseIndexManager` Class** (`src/core/InverseIndexManager.java`)
- Manages inverse permutation arrays for efficient range intersection
- Builds `inverseIndex[treeIndex][taxonId] = position` mappings
- Provides O(min(|A|, |B|)) intersection counting

#### 1.2 Inverse Index Weight Calculation Algorithm

**Core Algorithm:**
```java
// For intersection |A ∩ B| where A and B are ranges from different trees
int calculateRangeIntersection(int tree1, int start1, int end1, 
                              int tree2, int start2, int end2) {
    // Choose smaller range for iteration (complexity optimization)
    if ((end1 - start1) <= (end2 - start2)) {
        return countIntersection(tree1, start1, end1, tree2, start2, end2);
    } else {
        return countIntersection(tree2, start2, end2, tree1, start1, end1);
    }
}

private int countIntersection(int smallTree, int smallStart, int smallEnd,
                             int largeTree, int largeStart, int largeEnd) {
    int count = 0;
    int[] smallOrdering = geneTreeOrderings[smallTree];
    int[][] inverseIndex = this.inverseIndex;
    
    // Iterate over smaller range
    for (int pos = smallStart; pos < smallEnd; pos++) {
        int taxonId = smallOrdering[pos];
        int positionInLargeTree = inverseIndex[largeTree][taxonId];
        
        // Check if taxon falls within large tree's range
        if (positionInLargeTree >= largeStart && positionInLargeTree < largeEnd) {
            count++;
        }
    }
    return count;
}
```

**Complexity Analysis:**
- Per intersection: O(min(|A|, |B|)) instead of O(n)
- Total weight calculation: O(∑ min(|Ai|, |Bj|)) which is much smaller than O(n²k) under balanced tree assumptions

#### 1.3 New Weight Calculator Implementation

**`MemoryOptimizedWeightCalculator` Class** (`src/core/MemoryOptimizedWeightCalculator.java`)
- Uses range-based intersection counting via `InverseIndexManager`
- Supports CPU single-thread, CPU parallel, and GPU modes
- Works directly with `RangeBipartition` data from `MemoryEfficientBipartitionManager`
- No BitSet expansion required

### Phase 2: Memory-Optimized DP Algorithm

#### 2.1 Hash-Based Cluster Keys

**Problem:** Current DP uses `Map<BitSet, ...>` which requires O(n) memory per cluster key.

**Solution:** Replace with `Map<HashPair, ...>` using double-hashed cluster representations.

**`ClusterHashPair` Class** (`src/tree/ClusterHashPair.java`)
```java
public class ClusterHashPair {
    public final long sumHash;
    public final long xorHash;
    
    // Same structure as RangeBipartition.HashPair but for clusters
    // Represents hash of a contiguous range (i, start, end)
}
```

#### 2.2 Cluster Hash Computation

**`ClusterHashManager` Class** (`src/core/ClusterHashManager.java`)
```java
public class ClusterHashManager {
    private final Map<ClusterHashPair, Set<Integer>> hashToTaxonSets; // For verification
    private final HashUtils hashUtils;
    
    // Compute hash for a cluster defined by range
    public ClusterHashPair getClusterHash(int geneTreeIndex, int start, int end);
    
    // Compute hash for existing BitSet cluster (for compatibility)
    public ClusterHashPair getClusterHash(BitSet cluster);
    
    // Verify cluster equality with fallback to expensive verification
    public boolean clustersEqual(ClusterHashPair hash1, ClusterHashPair hash2);
    
    // Get left and right cluster hashes from a bipartition
    public ClusterHashPair getLeftClusterHash(RangeBipartition bip);
    public ClusterHashPair getRightClusterHash(RangeBipartition bip);
}
```

#### 2.3 DP State Management

**`MemoryOptimizedInferenceDP` Class** (`src/core/MemoryOptimizedInferenceDP.java`)
```java
public class MemoryOptimizedInferenceDP {
    private Map<ClusterHashPair, Double> dpMemo; // Hash-based memoization
    private Map<ClusterHashPair, STBipartition> dpChoice;
    private Map<ClusterHashPair, List<STBipartition>> clusterToSTBips; // Hash-based lookup
    private ClusterHashManager clusterHashManager;
    
    // Main DP function using cluster hashes
    private double dp(ClusterHashPair clusterHash);
    
    // Preprocessing: convert BitSet-based mappings to hash-based
    private void preprocessCandidatesWithHashes();
}
```

**Key DP Optimizations:**
1. **Hash-based memoization**: `dpMemo` uses `ClusterHashPair` keys instead of `BitSet`
2. **Hash-based candidate lookup**: `clusterToSTBips` maps cluster hashes to candidate bipartitions
3. **Lazy verification**: Trust double hash for equality, fallback to expensive checks only when needed

### Phase 3: GPU Optimization

#### 3.1 Compact GPU Data Transfer

**Current GPU Memory Usage:**
- Candidate BitSets: `numCandidates × bitsetSize × 8` bytes × 2 (left + right)
- Gene tree BitSets: `numGeneTreeBips × bitsetSize × 8` bytes × 2 (left + right)
- Total: O(n²k) memory transfer

**Optimized GPU Memory Usage:**
- Range representations: `numCandidates × 20` bytes (5 ints per range)
- Gene tree ranges: `numGeneTreeBips × 20` bytes  
- Inverse index arrays: `k × n × 4` bytes
- Permutation arrays: `k × n × 4` bytes
- Total: O(nk) memory transfer

#### 3.2 Updated CUDA Interface

**Modified `WeightCalcLib` Interface:**
```java
public interface WeightCalcLib extends Library {
    // New compact bipartition structure
    @Structure.FieldOrder({"geneTreeIndex", "leftStart", "leftEnd", "rightStart", "rightEnd"})
    public static class CompactBipartition extends Structure {
        public int geneTreeIndex;
        public int leftStart;
        public int leftEnd; 
        public int rightStart;
        public int rightEnd;
    }
    
    // New kernel launch method
    void launchCompactWeightCalculation(
        CompactBipartition[] candidates,
        CompactBipartition[] geneTreeBips,
        int[] frequencies,
        double[] weights,
        int[][] inverseIndex,      // [tree][taxon] = position
        int[][] geneTreeOrderings, // [tree][position] = taxon
        int numCandidates,
        int numGeneTreeBips,
        int numTrees,
        int numTaxa
    );
}
```

#### 3.3 CUDA Kernel Updates

**New GPU Intersection Algorithm:**
```cuda
__device__ int calculateCompactIntersection(
    CompactBipartition* bip1, CompactBipartition* bip2,
    int** inverseIndex, int** orderings) {
    
    // Choose smaller range for iteration
    int tree1 = bip1->geneTreeIndex;
    int tree2 = bip2->geneTreeIndex;
    
    // Calculate all four intersections: AA, AB, BA, BB
    int aa = rangeIntersection(tree1, bip1->leftStart, bip1->leftEnd,
                              tree2, bip2->leftStart, bip2->leftEnd,
                              inverseIndex, orderings);
    // ... similar for AB, BA, BB
    
    return calculateScore(aa, ab, ba, bb);
}
```

## Detailed Implementation Steps

### Step 1: Core Infrastructure Classes

#### 1.1 `Cluster` Class (`src/tree/Cluster.java`)
```java
package tree;

import utils.HashUtils;

public class Cluster {
    public final int geneTreeIndex;
    public final int startPos;
    public final int endPos;
    private ClusterHashPair cachedHash;
    
    public Cluster(int geneTreeIndex, int startPos, int endPos) {
        this.geneTreeIndex = geneTreeIndex;
        this.startPos = startPos;
        this.endPos = endPos;
    }
    
    public ClusterHashPair computeHash(HashUtils hashUtils, 
                                      int[][] geneTreeOrderings,
                                      long[][] prefixSums, 
                                      long[][] prefixXORs) {
        if (cachedHash == null) {
            cachedHash = hashUtils.computeClusterHash(
                geneTreeIndex, startPos, endPos,
                geneTreeOrderings, prefixSums, prefixXORs);
        }
        return cachedHash;
    }
    
    public int size() {
        return endPos - startPos;
    }
}
```

#### 1.2 `ClusterHashPair` Class (`src/tree/ClusterHashPair.java`)
```java
package tree;

public class ClusterHashPair {
    public final long sumHash;
    public final long xorHash;
    
    public ClusterHashPair(long sumHash, long xorHash) {
        this.sumHash = sumHash;
        this.xorHash = xorHash;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof ClusterHashPair)) return false;
        ClusterHashPair other = (ClusterHashPair) obj;
        return this.sumHash == other.sumHash && this.xorHash == other.xorHash;
    }
    
    @Override
    public int hashCode() {
        return (int)(sumHash ^ (sumHash >>> 32)) ^ (int)(xorHash ^ (xorHash >>> 32));
    }
    
    @Override
    public String toString() {
        return "ClusterHashPair(" + sumHash + ", " + xorHash + ")";
    }
}
```

#### 1.3 `HashUtils` Class (`src/utils/HashUtils.java`)
```java
package utils;

import tree.ClusterHashPair;

public class HashUtils {
    // Reuse RangeBipartition's hashing constants and methods
    private static final long MIX_CONST1 = 0x9e3779b97f4a7c15L;
    private static final long MIX_CONST2 = 0xc2b2ae3d27d4eb4fL;
    
    public static ClusterHashPair computeClusterHash(
            int geneTreeIndex, int start, int end,
            int[][] orderings, long[][] prefixSums, long[][] prefixXORs) {
        
        // Compute sum hash (similar to RangeBipartition.SumHashFunction)
        long sumHash = computeSumHash(geneTreeIndex, start, end, prefixSums);
        
        // Compute XOR hash (similar to RangeBipartition.XORHashFunction)  
        long xorHash = computeXORHash(geneTreeIndex, start, end, prefixXORs);
        
        return new ClusterHashPair(sumHash, xorHash);
    }
    
    private static long computeSumHash(int treeIndex, int start, int end, long[][] prefixSums) {
        if (prefixSums == null || treeIndex >= prefixSums.length) return 0L;
        long[] prefixSum = prefixSums[treeIndex];
        if (prefixSum == null) return 0L;
        
        long rangeSum = rangeSum(prefixSum, start, end);
        int size = end - start;
        
        // Apply mixing similar to RangeBipartition.SumHashFunction
        long combined = rangeSum * MIX_CONST1 ^ size;
        return Long.rotateLeft(combined, 27) ^ (combined >>> 33);
    }
    
    private static long computeXORHash(int treeIndex, int start, int end, long[][] prefixXORs) {
        if (prefixXORs == null || treeIndex >= prefixXORs.length) return 0L;
        long[] prefixXOR = prefixXORs[treeIndex];
        if (prefixXOR == null) return 0L;
        
        long rangeXOR = rangeXOR(prefixXOR, start, end);
        int size = end - start;
        
        // Apply mixing similar to RangeBipartition.XORHashFunction
        return xorClusterMix(rangeXOR, size, MIX_CONST2);
    }
    
    // Helper methods (similar to RangeBipartition implementations)
    private static long rangeSum(long[] prefix, int start, int end) {
        if (end == 0) return 0;
        long startSum = (start > 0) ? prefix[start - 1] : 0;
        return prefix[end - 1] - startSum;
    }
    
    private static long rangeXOR(long[] prefixXOR, int start, int end) {
        if (prefixXOR == null || end <= 0) return 0L;
        int e = Math.min(end - 1, prefixXOR.length - 1);
        long endXOR = prefixXOR[e];
        long startXOR = (start > 0) ? prefixXOR[Math.min(start - 1, prefixXOR.length - 1)] : 0L;
        return endXOR ^ startXOR;
    }
    
    private static long xorClusterMix(long xorVal, int size, long salt) {
        long v = xorVal ^ (((long) size) << 32) ^ salt;
        return splitMix64(v);
    }
    
    private static long splitMix64(long z) {
        z += 0x9e3779b97f4a7c15L;
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }
}
```

### Step 2: Inverse Index Manager

#### 2.1 `InverseIndexManager` Class (`src/core/InverseIndexManager.java`)
```java
package core;

import java.util.List;
import tree.Tree;
import tree.TreeNode;

public class InverseIndexManager {
    private final int[][] inverseIndex; // [treeIndex][taxonId] = position
    private final int[][] geneTreeOrderings; // [treeIndex][position] = taxonId
    private final int numTrees;
    private final int numTaxa;
    
    public InverseIndexManager(List<Tree> geneTrees, int realTaxaCount) {
        this.numTrees = geneTrees.size();
        this.numTaxa = realTaxaCount;
        this.inverseIndex = new int[numTrees][numTaxa];
        this.geneTreeOrderings = new int[numTrees][];
        
        buildInverseIndex(geneTrees);
    }
    
    private void buildInverseIndex(List<Tree> geneTrees) {
        System.out.println("Building inverse index for " + numTrees + " trees...");
        
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            Tree tree = geneTrees.get(treeIdx);
            
            // Get left-to-right ordering (same as MemoryEfficientBipartitionManager)
            List<Integer> ordering = new ArrayList<>();
            collectLeavesInOrder(tree.root, ordering);
            
            geneTreeOrderings[treeIdx] = ordering.stream().mapToInt(Integer::intValue).toArray();
            
            // Build inverse mapping: taxonId -> position
            for (int pos = 0; pos < geneTreeOrderings[treeIdx].length; pos++) {
                int taxonId = geneTreeOrderings[treeIdx][pos];
                inverseIndex[treeIdx][taxonId] = pos;
            }
        }
        
        System.out.println("Inverse index built successfully");
    }
    
    private void collectLeavesInOrder(TreeNode node, List<Integer> ordering) {
        if (node.isLeaf()) {
            ordering.add(node.taxon.id);
            return;
        }
        
        if (node.childs != null && node.childs.size() >= 2) {
            collectLeavesInOrder(node.childs.get(0), ordering);
            collectLeavesInOrder(node.childs.get(1), ordering);
        }
        
        for (int i = 2; i < (node.childs != null ? node.childs.size() : 0); i++) {
            collectLeavesInOrder(node.childs.get(i), ordering);
        }
    }
    
    /**
     * Calculate intersection size between two ranges using inverse index.
     * Complexity: O(min(|range1|, |range2|))
     */
    public int getRangeIntersectionSize(int tree1, int start1, int end1, 
                                       int tree2, int start2, int end2) {
        // Choose smaller range for iteration
        int size1 = end1 - start1;
        int size2 = end2 - start2;
        
        if (size1 <= size2) {
            return countIntersection(tree1, start1, end1, tree2, start2, end2);
        } else {
            return countIntersection(tree2, start2, end2, tree1, start1, end1);
        }
    }
    
    private int countIntersection(int smallTree, int smallStart, int smallEnd,
                                 int largeTree, int largeStart, int largeEnd) {
        int count = 0;
        int[] smallOrdering = geneTreeOrderings[smallTree];
        
        for (int pos = smallStart; pos < smallEnd && pos < smallOrdering.length; pos++) {
            int taxonId = smallOrdering[pos];
            int positionInLargeTree = inverseIndex[largeTree][taxonId];
            
            if (positionInLargeTree >= largeStart && positionInLargeTree < largeEnd) {
                count++;
            }
        }
        
        return count;
    }
    
    // Getters for GPU implementation
    public int[][] getInverseIndex() { return inverseIndex; }
    public int[][] getGeneTreeOrderings() { return geneTreeOrderings; }
}
```

### Step 3: Memory-Optimized Weight Calculator

#### 3.1 `MemoryOptimizedWeightCalculator` Class (`src/core/MemoryOptimizedWeightCalculator.java`)
```java
package core;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;

import preprocessing.GeneTrees;
import tree.STBipartition;
import tree.RangeBipartition;
import tree.MemoryEfficientBipartitionManager;
import utils.Config;
import utils.Threading;

public class MemoryOptimizedWeightCalculator {
    
    private final GeneTrees geneTrees;
    private final InverseIndexManager inverseIndexManager;
    private final Map<Object, List<RangeBipartition>> hashToBipartitions;
    
    public MemoryOptimizedWeightCalculator(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
        this.inverseIndexManager = new InverseIndexManager(geneTrees.geneTrees, geneTrees.realTaxaCount);
        
        // Get range bipartitions from existing MemoryEfficientBipartitionManager
        MemoryEfficientBipartitionManager manager = 
            new MemoryEfficientBipartitionManager(geneTrees.geneTrees, geneTrees.realTaxaCount);
        manager.processGeneTreesParallel(); // This populates the hash mappings
        this.hashToBipartitions = manager.getHashToBipartitions(); // Need to expose this
    }
    
    public Map<STBipartition, Double> calculateWeights(List<STBipartition> candidates) {
        System.out.println("==== MEMORY-OPTIMIZED WEIGHT CALCULATION STARTED ====");
        
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
    
    private Map<STBipartition, Double> calculateWeightsSingleThread(List<STBipartition> candidates) {
        Map<STBipartition, Double> weights = new HashMap<>();
        
        // Convert candidates to range representations for efficient processing
        Map<STBipartition, RangeBipartition> candidateToRange = buildCandidateRangeMapping(candidates);
        
        for (STBipartition candidate : candidates) {
            double totalScore = 0.0;
            RangeBipartition candidateRange = candidateToRange.get(candidate);
            
            // Iterate through gene tree bipartitions using range representations
            for (Map.Entry<Object, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
                List<RangeBipartition> geneTreeRanges = entry.getValue();
                int frequency = geneTreeRanges.size(); // All ranges in same hash group
                
                if (!geneTreeRanges.isEmpty()) {
                    RangeBipartition geneTreeRange = geneTreeRanges.get(0); // Representative
                    double score = calculateRangeScore(candidateRange, geneTreeRange);
                    totalScore += score * frequency;
                }
            }
            
            weights.put(candidate, totalScore);
        }
        
        return weights;
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
    
    private Map<STBipartition, RangeBipartition> buildCandidateRangeMapping(List<STBipartition> candidates) {
        // This would need to convert STBipartitions back to range representations
        // Implementation depends on how we want to handle this conversion
        // For now, placeholder - would need to be implemented based on specific needs
        return new HashMap<>();
    }
}
```

### Step 4: Memory-Optimized DP Implementation

#### 4.1 `ClusterHashManager` Class (`src/core/ClusterHashManager.java`)
```java
package core;

import java.util.*;
import tree.ClusterHashPair;
import tree.RangeBipartition;
import utils.BitSet;
import utils.HashUtils;

public class ClusterHashManager {
    private final Map<ClusterHashPair, Set<Integer>> hashToTaxonSets;
    private final int[][] geneTreeOrderings;
    private final long[][] prefixSums;
    private final long[][] prefixXORs;
    private final int realTaxaCount;
    
    public ClusterHashManager(int[][] geneTreeOrderings, long[][] prefixSums, 
                             long[][] prefixXORs, int realTaxaCount) {
        this.hashToTaxonSets = new HashMap<>();
        this.geneTreeOrderings = geneTreeOrderings;
        this.prefixSums = prefixSums;
        this.prefixXORs = prefixXORs;
        this.realTaxaCount = realTaxaCount;
    }
    
    /**
     * Compute hash for a cluster defined by range in a gene tree.
     */
    public ClusterHashPair getClusterHash(int geneTreeIndex, int start, int end) {
        return HashUtils.computeClusterHash(geneTreeIndex, start, end,
                                           geneTreeOrderings, prefixSums, prefixXORs);
    }
    
    /**
     * Compute hash for existing BitSet cluster (for compatibility).
     */
    public ClusterHashPair getClusterHash(BitSet cluster) {
        // Convert BitSet to a canonical range representation for hashing
        // This is a fallback method for compatibility with existing code
        
        // Find a representative gene tree and range that matches this cluster
        // This is expensive but needed for backward compatibility
        Set<Integer> taxonSet = new HashSet<>();
        for (int i = cluster.nextSetBit(0); i >= 0; i = cluster.nextSetBit(i + 1)) {
            taxonSet.add(i);
        }
        
        // Try to find existing hash for this taxon set
        for (Map.Entry<ClusterHashPair, Set<Integer>> entry : hashToTaxonSets.entrySet()) {
            if (entry.getValue().equals(taxonSet)) {
                return entry.getKey();
            }
        }
        
        // If not found, compute new hash (expensive path)
        // For now, use a deterministic hash based on taxon set content
        return computeFallbackHash(taxonSet);
    }
    
    /**
     * Get left cluster hash from a bipartition.
     */
    public ClusterHashPair getLeftClusterHash(RangeBipartition bip) {
        ClusterHashPair hash = getClusterHash(bip.geneTreeIndex, bip.leftStart, bip.leftEnd);
        
        // Cache the taxon set for verification if needed
        Set<Integer> taxonSet = getRangeTaxonSet(bip.geneTreeIndex, bip.leftStart, bip.leftEnd);
        hashToTaxonSets.put(hash, taxonSet);
        
        return hash;
    }
    
    /**
     * Get right cluster hash from a bipartition.
     */
    public ClusterHashPair getRightClusterHash(RangeBipartition bip) {
        ClusterHashPair hash = getClusterHash(bip.geneTreeIndex, bip.rightStart, bip.rightEnd);
        
        // Cache the taxon set for verification if needed
        Set<Integer> taxonSet = getRangeTaxonSet(bip.geneTreeIndex, bip.rightStart, bip.rightEnd);
        hashToTaxonSets.put(hash, taxonSet);
        
        return hash;
    }
    
    /**
     * Verify cluster equality with fallback to expensive verification.
     */
    public boolean clustersEqual(ClusterHashPair hash1, ClusterHashPair hash2) {
        if (hash1.equals(hash2)) {
            return true; // Trust the double hash
        }
        
        // Hash collision - need expensive verification
        Set<Integer> set1 = hashToTaxonSets.get(hash1);
        Set<Integer> set2 = hashToTaxonSets.get(hash2);
        
        if (set1 != null && set2 != null) {
            return set1.equals(set2);
        }
        
        return false; // Cannot verify - assume different
    }
    
    private Set<Integer> getRangeTaxonSet(int geneTreeIndex, int start, int end) {
        Set<Integer> taxonSet = new HashSet<>();
        int[] ordering = geneTreeOrderings[geneTreeIndex];
        
        for (int i = start; i < end && i < ordering.length; i++) {
            taxonSet.add(ordering[i]);
        }
        
        return taxonSet;
    }
    
    private ClusterHashPair computeFallbackHash(Set<Integer> taxonSet) {
        // Compute a deterministic hash for the taxon set
        // This is a fallback when we can't find a range representation
        
        long sumHash = 0;
        long xorHash = 0;
        
        for (int taxonId : taxonSet) {
            long hashedTaxon = hashSingleTaxon(taxonId);
            sumHash += hashedTaxon;
            xorHash ^= hashedTaxon;
        }
        
        // Apply mixing
        sumHash = Long.rotateLeft(sumHash * 0x9e3779b97f4a7c15L, 27);
        xorHash = splitMix64(xorHash ^ (taxonSet.size() << 16));
        
        return new ClusterHashPair(sumHash, xorHash);
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
    
    private static long splitMix64(long z) {
        z += 0x9e3779b97f4a7c15L;
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }
}
```

#### 4.2 `MemoryOptimizedInferenceDP` Class (`src/core/MemoryOptimizedInferenceDP.java`)
```java
package core;

import java.util.*;
import preprocessing.GeneTrees;
import tree.*;
import utils.BitSet;

public class MemoryOptimizedInferenceDP {
    
    private final GeneTrees geneTrees;
    private final List<STBipartition> candidateSTBips;
    private final Map<ClusterHashPair, List<STBipartition>> clusterHashToSTBips;
    private final Map<STBipartition, Double> stbipWeights;
    private final Map<ClusterHashPair, Double> dpMemo;
    private final Map<ClusterHashPair, STBipartition> dpChoice;
    private final ClusterHashManager clusterHashManager;
    
    public MemoryOptimizedInferenceDP(GeneTrees geneTrees, List<STBipartition> candidateSTBips) {
        this.geneTrees = geneTrees;
        this.candidateSTBips = candidateSTBips;
        this.clusterHashToSTBips = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
        
        // Initialize cluster hash manager with gene tree data
        // This would need access to the ordering and prefix arrays from MemoryEfficientBipartitionManager
        this.clusterHashManager = initializeClusterHashManager();
        
        preprocessCandidatesWithHashes();
        calculateWeights();
    }
    
    private ClusterHashManager initializeClusterHashManager() {
        // Get the data structures from MemoryEfficientBipartitionManager
        // This would require exposing these arrays from the manager
        MemoryEfficientBipartitionManager manager = 
            new MemoryEfficientBipartitionManager(geneTrees.geneTrees, geneTrees.realTaxaCount);
        
        // Would need to expose these from MemoryEfficientBipartitionManager:
        int[][] geneTreeOrderings = manager.getGeneTreeOrderings();
        long[][] prefixSums = manager.getPrefixSums();
        long[][] prefixXORs = manager.getPrefixXORs();
        
        return new ClusterHashManager(geneTreeOrderings, prefixSums, prefixXORs, geneTrees.realTaxaCount);
    }
    
    /**
     * Preprocess candidates to create hash-based cluster mappings.
     */
    private void preprocessCandidatesWithHashes() {
        System.out.println("Preprocessing candidates with cluster hashes...");
        
        for (STBipartition stbip : candidateSTBips) {
            // Compute union cluster hash
            BitSet union = (BitSet) stbip.cluster1.clone();
            union.or(stbip.cluster2);
            
            ClusterHashPair unionHash = clusterHashManager.getClusterHash(union);
            clusterHashToSTBips.computeIfAbsent(unionHash, k -> new ArrayList<>()).add(stbip);
        }
        
        System.out.println("Created hash mappings for " + clusterHashToSTBips.size() + " unique clusters");
    }
    
    private void calculateWeights() {
        MemoryOptimizedWeightCalculator calculator = new MemoryOptimizedWeightCalculator(geneTrees);
        stbipWeights = calculator.calculateWeights(candidateSTBips);
    }
    
    public double solve() {
        BitSet allTaxa = new BitSet(geneTrees.realTaxaCount);
        for (int i = 0; i < geneTrees.realTaxaCount; i++) {
            allTaxa.set(i);
        }
        
        ClusterHashPair allTaxaHash = clusterHashManager.getClusterHash(allTaxa);
        return dp(allTaxaHash);
    }
    
    /**
     * Main DP function using cluster hashes instead of BitSets.
     */
    private double dp(ClusterHashPair clusterHash) {
        if (dpMemo.containsKey(clusterHash)) {
            return dpMemo.get(clusterHash);
        }
        
        // Get cluster size - this requires looking up the actual cluster
        // For now, we'll need a way to get cluster size from hash
        int taxaCount = getClusterSize(clusterHash);
        
        if (taxaCount <= 2) {
            dpMemo.put(clusterHash, 0.0);
            return 0.0;
        }
        
        double maxScore = Double.NEGATIVE_INFINITY;
        STBipartition bestChoice = null;
        
        List<STBipartition> candidates = clusterHashToSTBips.get(clusterHash);
        if (candidates != null) {
            for (STBipartition stbip : candidates) {
                if (isValidPartition(stbip, clusterHash)) {
                    // Get left and right cluster hashes
                    ClusterHashPair leftHash = clusterHashManager.getClusterHash(stbip.cluster1);
                    ClusterHashPair rightHash = clusterHashManager.getClusterHash(stbip.cluster2);
                    
                    double leftScore = dp(leftHash);
                    double rightScore = dp(rightHash);
                    double stbipScore = stbipWeights.getOrDefault(stbip, 0.0);
                    
                    double totalScore = leftScore + rightScore + stbipScore;
                    
                    if (totalScore > maxScore) {
                        maxScore = totalScore;
                        bestChoice = stbip;
                    }
                }
            }
        }
        
        if (bestChoice == null) {
            maxScore = 0.0;
        }
        
        dpMemo.put(clusterHash, maxScore);
        dpChoice.put(clusterHash, bestChoice);
        return maxScore;
    }
    
    private boolean isValidPartition(STBipartition stbip, ClusterHashPair clusterHash) {
        // Verify that stbip.cluster1 ∪ stbip.cluster2 == cluster represented by clusterHash
        BitSet union = (BitSet) stbip.cluster1.clone();
        union.or(stbip.cluster2);
        
        ClusterHashPair unionHash = clusterHashManager.getClusterHash(union);
        return clusterHashManager.clustersEqual(unionHash, clusterHash);
    }
    
    private int getClusterSize(ClusterHashPair clusterHash) {
        // This is a challenge - we need to get cluster size from hash
        // Options:
        // 1. Cache size information in ClusterHashManager
        // 2. Store size as part of ClusterHashPair
        // 3. Look up actual cluster from hash (expensive)
        
        // For now, placeholder - would need to be implemented
        // Could cache this information during preprocessing
        return 0; // TODO: Implement proper cluster size lookup
    }
    
    // Rest of the methods (reconstructTree, buildTreeNode, etc.) would be similar
    // but work with cluster hashes instead of BitSets where possible
}
```

### Step 5: GPU Optimization

#### 5.1 Updated CUDA Interface (`src/core/WeightCalculator.java`)
```java
// Add to existing WeightCalcLib interface
@Structure.FieldOrder({"geneTreeIndex", "leftStart", "leftEnd", "rightStart", "rightEnd"})
public static class CompactBipartition extends Structure {
    public int geneTreeIndex;
    public int leftStart;
    public int leftEnd;
    public int rightStart;
    public int rightEnd;
    
    public CompactBipartition() {
        super();
    }
}

void launchCompactWeightCalculation(
    CompactBipartition[] candidates,
    CompactBipartition[] geneTreeBips,
    int[] frequencies,
    double[] weights,
    Pointer inverseIndexPtr,    // Flattened [tree*numTaxa + taxon] = position
    Pointer orderingPtr,        // Flattened [tree*numTaxa + position] = taxon
    int numCandidates,
    int numGeneTreeBips,
    int numTrees,
    int numTaxa
);
```

#### 5.2 Compact GPU Weight Calculator
```java
private Map<STBipartition, Double> calculateWeightsCompactGPU(List<STBipartition> candidates) {
    // Convert STBipartitions to CompactBipartitions
    // Transfer compact data structures instead of full BitSets
    // Memory usage: O(nk) instead of O(n²k)
    
    System.out.println("==== COMPACT GPU WEIGHT CALCULATION ====");
    
    // Much smaller memory allocations
    CompactBipartition[] compactCandidates = convertToCompactBipartitions(candidates);
    CompactBipartition[] compactGeneTreeBips = convertGeneTreeBipartitions();
    
    // Transfer inverse index and ordering arrays
    Memory inverseIndexMem = flattenInverseIndex();
    Memory orderingMem = flattenOrderings();
    
    // Launch compact kernel
    WeightCalcLib.INSTANCE.launchCompactWeightCalculation(
        compactCandidates, compactGeneTreeBips, frequencies, weights,
        inverseIndexMem, orderingMem, numCandidates, numGeneTreeBips, 
        numTrees, numTaxa);
    
    return convertResults(weights, candidates);
}
```

## Key Design Principles

### 1. **Isolation and Compatibility**
- New classes don't modify existing `RangeBipartition` or `MemoryEfficientBipartitionManager`
- Maintain existing interfaces and add new optimized paths
- Backward compatibility with existing computation modes

### 2. **Hash Consistency**
- **Double Hashing**: All cluster hashes use `ClusterHashPair` (sum hash, xor hash) for collision resistance
- **Same Scheme**: Reuse the proven double-hashing approach from `RangeBipartition`
- **Isolation**: Cluster hashing is separate from bipartition hashing but uses same underlying algorithms

### 3. **Memory Optimization Strategy**
- **Range-Based Processing**: Work with compact integer tuples throughout the pipeline
- **Lazy Expansion**: Only convert to BitSets when absolutely necessary
- **Hash-Based Lookups**: Use double hashes as keys instead of full data structures

### 4. **Verification and Trust**
- **Trust Double Hash**: Assume double hash collisions are extremely rare
- **Fallback Verification**: Expensive equality checks available when needed
- **Configurable**: Can enable/disable expensive verification for debugging

### 5. **Modularity**
- **Independent Optimization**: Weight calculation, DP, and GPU can be optimized separately
- **Pluggable Components**: Easy to switch between old and new implementations
- **Incremental Deployment**: Can optimize one component at a time

## Expected Benefits

### 1. **Memory Reduction**
- **Weight Calculation**: O(nk) instead of O(n²k) memory usage
- **DP Algorithm**: O(number of unique clusters) instead of O(n × number of clusters)  
- **GPU Transfer**: O(nk) data transfer instead of O(n²k)
- **Overall**: Significant reduction in peak memory usage

### 2. **Performance Improvements**
- **Faster Intersections**: O(min(|A|, |B|)) per intersection instead of O(n)
- **Reduced Allocations**: Less memory allocation and garbage collection pressure
- **Better Cache Locality**: Compact data structures improve cache performance
- **GPU Efficiency**: Better GPU utilization with compact data transfers

### 3. **Scalability**
- **Larger Datasets**: Handle datasets that previously caused out-of-memory errors
- **Memory-Constrained Systems**: Better performance on systems with limited RAM
- **Future-Proof**: Architecture supports even larger datasets as they become available

## Implementation Timeline

### Phase 1: Core Infrastructure (Week 1-2)
- Implement `Cluster`, `ClusterHashPair`, `HashUtils` classes
- Create `InverseIndexManager` with range intersection algorithms
- Unit tests for core functionality

### Phase 2: Weight Calculator Optimization (Week 3-4)  
- Implement `MemoryOptimizedWeightCalculator`
- CPU single-thread and parallel modes
- Integration with existing `MemoryEfficientBipartitionManager`
- Performance benchmarking vs. original implementation

### Phase 3: DP Optimization (Week 5-6)
- Implement `ClusterHashManager` and `MemoryOptimizedInferenceDP`
- Hash-based cluster operations and memoization
- Validation against original DP results

### Phase 4: GPU Optimization (Week 7-8)
- Update CUDA interface and kernel for compact data
- Implement compact GPU weight calculator
- Performance testing and optimization

### Phase 5: Integration and Testing (Week 9-10)
- Full system integration
- Comprehensive testing on various datasets
- Performance benchmarking and memory profiling
- Documentation and code review

## Risk Mitigation

### 1. **Hash Collision Handling**
- **Double Hashing**: Use both sum and XOR hashes to minimize collision probability
- **Verification Framework**: Expensive equality checks available as fallback
- **Monitoring**: Log and track hash collision rates during testing

### 2. **Correctness Validation**
- **Cross-Validation**: Compare results between old and new implementations
- **Unit Tests**: Comprehensive test coverage for all new components
- **Integration Tests**: End-to-end validation on known datasets

### 3. **Performance Regression**
- **Benchmarking**: Measure performance at each implementation phase
- **Profiling**: Memory and CPU profiling to identify bottlenecks
- **Fallback Options**: Keep original implementations available if needed

### 4. **Memory Management**
- **Careful Allocation**: Monitor memory usage during development
- **Garbage Collection**: Minimize object creation in hot paths
- **Resource Cleanup**: Proper cleanup of native memory in GPU code

This comprehensive plan provides a roadmap for implementing memory-optimized weight calculation and DP algorithms while maintaining the robustness and correctness of the existing system. The double-hashing approach for clusters ensures collision resistance while providing the memory benefits of compact hash-based representations.
