# ASTRAL Bipartition Set Expansion Analysis

## Overview

ASTRAL uses a sophisticated multi-stage approach to expand the set of bipartitions (called "set X") beyond just the bipartitions found in input gene trees. This expansion is crucial for finding optimal species trees, especially when gene trees are incomplete or have missing data.

## ASTRAL's Bipartition Expansion Strategy

### Stage 1: Basic Gene Tree Bipartitions
- **Location**: `addBipartitionsFromSignleIndTreesToX()` method
- **Purpose**: Extract all bipartitions from input gene trees
- **Process**: 
  - Traverse each gene tree in post-order
  - For each internal node, create a bipartition representing the split
  - Add bipartition to set X if it's valid (not trivial)

### Stage 2: Distance-Based Bipartitions  
- **Location**: `addExtraBipartitionByDistance()` method (line 1029)
- **Purpose**: Add bipartitions from distance-based trees (UPGMA/NJ)
- **Process**:
  ```java
  // Build distance matrix from gene trees
  calculateDistances();
  
  // Infer tree from distance matrix and extract bipartitions
  for (BitSet bs : speciesMatrix.inferTreeBitsets()) {
      STITreeCluster g = GlobalMaps.taxonNameMap.getSpeciesIdMapper()
              .getGeneClusterForSTCluster(bs);
      this.addCompletedSpeciesFixedBipartionToX(g, g.complementaryCluster());
  }
  
  // In SLOW mode, also add quadratic bipartitions
  if (SLOW) {
      for (BitSet bs : speciesMatrix.getQuadraticBitsets()) {
          // Add more bipartitions from quartet-based methods
      }
  }
  ```

### Stage 3: Greedy Consensus Bipartitions
- **Location**: `addExtraBipartitionByHeuristics()` method (line 1061)
- **Purpose**: Generate consensus trees and extract their bipartitions
- **Key Innovation**: Creates multiple greedy consensus trees with different thresholds

#### Greedy Consensus Process:
```java
// Create greedy consensus trees with multiple thresholds
allGreedies = Utils.greedyConsensus(contractedTrees,
        this.GREEDY_ADDITION_THRESHOLDS, true, 1, tid, true);

// GREEDY_ADDITION_THRESHOLDS = {0, 1/100, 1/50, 1/20, 1/10, 1/5, 1/3}
// This creates 7 different consensus trees with increasing strictness
```

#### Polytomy Resolution:
- **Problem**: Greedy consensus trees often have polytomies (nodes with >2 children)
- **Solution**: Resolve polytomies using multiple strategies:

1. **Distance-based resolution** (line 1234):
   ```java
   // Resolve polytomy using distance matrix
   WQDataCollection.this.addSubSampledBitSetToX(
       WQDataCollection.this.speciesMatrix.resolvePolytomy(
           Arrays.asList(childbs), true), tid);
   ```

2. **Sampling-based resolution** (line 1250):
   ```java
   // Sample taxa around polytomy and resolve using greedy consensus
   if (sampleAndResolve(childbs, contractedTrees, quadratic, tid, true, false)) {
       k += GREEDY_ADDITION_IMPROVEMENT_REWARD;
   }
   ```

### Stage 4: Multi-Round Sampling
- **Location**: Lines 667-811 in `formSetX()`
- **Purpose**: Handle multi-individual datasets through sampling
- **Process**:
  1. Create K=100 random samples per gene tree
  2. For each sample, create contracted trees
  3. Build greedy consensus from subsamples
  4. Extract bipartitions from all consensus trees

## Key Algorithmic Insights

### 1. Threshold-Based Consensus
ASTRAL creates multiple consensus trees using different support thresholds:
- Threshold 0: Include all bipartitions (most permissive)
- Threshold 1/3: Include bipartitions in ≥33% of gene trees (most restrictive)
- This captures bipartitions with varying levels of support

### 2. Polytomy Resolution Strategies
When consensus trees have polytomies, ASTRAL uses:
- **Distance methods**: UPGMA/NJ to resolve based on evolutionary distances
- **Sampling methods**: Random sampling around polytomy branches
- **Iterative improvement**: Reward successful resolutions with more attempts

### 3. Adaptive Polytomy Size Limits
```java
// Limit polytomy resolution based on computational complexity
int N = GREEDY_ADDITION_MAX_POLYTOMY_MIN + 
        speciesCount * GREEDY_ADDITION_MAX_POLYTOMY_MULT;
// Only resolve polytomies up to a certain size to maintain efficiency
```

## Application to STELAR-MP

### Current STELAR Approach
Your current STELAR implementation:
1. Extracts bipartitions directly from gene trees
2. Runs weight calculation on these bipartitions
3. Uses DP to find optimal species tree

### Potential Enhancements

#### 1. **Distance-Based Expansion** (Easiest to implement)
```java
// In your GeneTrees preprocessing
public void addDistanceBasedBipartitions() {
    // Build distance matrix from gene trees
    DistanceMatrix distMatrix = buildDistanceMatrix();
    
    // Infer UPGMA/NJ tree
    Tree distanceTree = distMatrix.buildUPGMATree();
    
    // Extract bipartitions from distance tree
    List<STBipartition> distanceBips = extractBipartitions(distanceTree);
    
    // Add to candidate set
    candidateSTBips.addAll(distanceBips);
}
```

#### 2. **Greedy Consensus Expansion** (Moderate complexity)
```java
// Create consensus trees from gene trees
public void addConsensusBasedBipartitions() {
    double[] thresholds = {0.0, 0.1, 0.2, 0.33, 0.5};
    
    for (double threshold : thresholds) {
        Tree consensus = buildGreedyConsensus(geneTrees, threshold);
        List<STBipartition> consensusBips = extractBipartitions(consensus);
        candidateSTBips.addAll(consensusBips);
    }
}
```

#### 3. **Polytomy Resolution** (Advanced)
```java
// Resolve polytomies in consensus trees
public void resolvePolytomies(Tree consensusTree) {
    for (TreeNode node : consensusTree.getNodes()) {
        if (node.getChildCount() > 2) {
            // Method 1: Distance-based resolution
            List<STBipartition> resolved = resolveByDistance(node);
            
            // Method 2: Sampling-based resolution  
            List<STBipartition> sampled = resolveBySampling(node);
            
            candidateSTBips.addAll(resolved);
            candidateSTBips.addAll(sampled);
        }
    }
}
```

### Implementation Strategy for STELAR

#### Phase 1: Distance-Based Enhancement
1. **Modify `GeneTrees` class**:
   - Add distance matrix calculation
   - Implement UPGMA/NJ tree construction
   - Extract bipartitions from distance trees

2. **Integration point**: 
   - Add distance-based bipartitions in `GeneTrees.preprocessTrees()`
   - Before calling `InferenceDP`

#### Phase 2: Consensus-Based Enhancement  
1. **Add consensus tree builder**:
   - Implement greedy consensus algorithm
   - Support multiple thresholds
   - Handle bipartition frequency counting

2. **Integration point**:
   - Add consensus-based bipartitions after distance-based ones
   - May significantly increase candidate set size

#### Phase 3: Advanced Polytomy Resolution
1. **Implement polytomy detection and resolution**
2. **Add sampling strategies for large polytomies**
3. **Implement adaptive size limits for computational efficiency**

### Expected Benefits

1. **Improved Accuracy**: More comprehensive bipartition set should lead to better species trees
2. **Robustness**: Distance and consensus methods can recover signal missed by individual gene trees
3. **Handling Missing Data**: Consensus methods are more robust to incomplete gene trees

### Computational Considerations

1. **Increased Candidate Set Size**: May grow from O(n×g) to O(n×g×k) where k is expansion factor
2. **Weight Calculation Impact**: Your GPU implementation should handle larger candidate sets well
3. **Memory Usage**: Need to monitor memory consumption with expanded bipartition sets

### Recommended Implementation Order

1. **Start with distance-based expansion** (lowest risk, moderate benefit)
2. **Add simple consensus with 2-3 thresholds** (moderate risk, high benefit)  
3. **Consider polytomy resolution if needed** (high complexity, incremental benefit)

This approach would make STELAR more similar to ASTRAL-II's comprehensive bipartition search while maintaining your efficient DP-based inference and GPU acceleration advantages.

## Summary

ASTRAL's bipartition expansion is a sophisticated multi-stage process that goes far beyond just extracting bipartitions from gene trees. The key insight is that the optimal species tree may contain bipartitions that don't appear in any individual gene tree but emerge from:

1. **Distance-based inference** - Capturing overall evolutionary signal
2. **Consensus methods** - Finding bipartitions supported across multiple gene trees
3. **Polytomy resolution** - Systematically exploring alternative resolutions of uncertain relationships

For STELAR-MP, implementing even a subset of these techniques (starting with distance-based expansion) could significantly improve accuracy while leveraging your existing GPU-accelerated weight calculation infrastructure.

The commented ASTRAL code now clearly shows the algorithmic flow and can serve as a reference for implementing similar techniques in your STELAR codebase.
