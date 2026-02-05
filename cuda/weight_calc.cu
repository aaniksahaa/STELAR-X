#include <stdio.h>
#include <cuda_runtime.h>

// Legacy structure to represent a bipartition using BitSets
struct Bipartition {
    unsigned long long* cluster1;  // First cluster bits
    unsigned long long* cluster2;  // Second cluster bits
    int bitsetSize;               // Size of bitset in 64-bit words
};

// New compact structure for memory-optimized bipartitions
struct CompactBipartition {
    int geneTreeIndex;    // Which gene tree this bipartition belongs to
    int leftStart;        // Start index of left subtree range (inclusive)
    int leftEnd;          // End index of left subtree range (exclusive)
    int rightStart;       // Start index of right subtree range (inclusive)
    int rightEnd;         // End index of right subtree range (exclusive)
};

// Mixed bipartition structure for cross-tree recombination
// Left and right sides can come from different gene trees
struct MixedCompactBipartition {
    int leftTreeIndex;    // Gene tree index for left side
    int leftStart;        // Start index of left subtree range (inclusive)
    int leftEnd;          // End index of left subtree range (exclusive)
    int rightTreeIndex;   // Gene tree index for right side (may differ from left)
    int rightStart;       // Start index of right subtree range (inclusive)
    int rightEnd;         // End index of right subtree range (exclusive)
};

// Structure for quartet scoring additional data
// Contains per-tree presence counts for efficient restricted size calculation
struct QuartetScoringData {
    int* presentCount;     // presentCount[treeIndex] = |Sg| (number of taxa in gene tree g)
    int* missingTaxaFlat;  // Flattened missing taxa list for all trees
    int* missingOffsets;   // missingOffsets[treeIndex] = start offset in missingTaxaFlat
    int* missingCounts;    // missingCounts[treeIndex] = number of missing taxa in tree g
};

// ASTRAL-style F function: F(a,b,c) = a * b * c * (a + b + c - 3)
// Note: No division by 2 - matches ASTRAL's actual implementation
__device__ double quartetF(int a, int b, int c) {
    if (a < 0 || b < 0 || c < 0) {
        return 0.0;
    }
    long sum = (long)a + b + c;
    if (sum < 3) {
        return 0.0;
    }
    // Use long arithmetic to avoid overflow - NO division by 2
    return (double)((long)a * b * c * (sum - 3));
}

// Compute QI score from a 3x3 intersection grid
// QI = sum of 6 terms: F(n00,n11,n22) + F(n00,n12,n21) + F(n01,n10,n22) + F(n01,n12,n20) + F(n02,n10,n21) + F(n02,n11,n20)
__device__ double computeQuartetQI(int* grid) {
    // grid is flattened [row * 3 + col]
    double qi = quartetF(grid[0], grid[4], grid[8])   // F(n00, n11, n22)
              + quartetF(grid[0], grid[5], grid[7])   // F(n00, n12, n21)
              + quartetF(grid[1], grid[3], grid[8])   // F(n01, n10, n22)
              + quartetF(grid[1], grid[5], grid[6])   // F(n01, n12, n20)
              + quartetF(grid[2], grid[3], grid[7])   // F(n02, n10, n21)
              + quartetF(grid[2], grid[4], grid[6]);  // F(n02, n11, n20)
    
    return qi;  // No division - matches ASTRAL's formula
}

// Count missing taxa from a candidate range that are not in the target gene tree
// Uses the missing taxa list for efficiency when missing count is small
__device__ int countMissingInRange(
    int sourceTree, int rangeStart, int rangeEnd,
    int targetTree,
    int* inverseIndex,
    int* orderings,
    int* missingTaxaFlat,
    int* missingOffsets,
    int* missingCounts,
    int numTaxa
) {
    int count = 0;
    int missingStart = missingOffsets[targetTree];
    int missingEnd = missingStart + missingCounts[targetTree];
    
    // Iterate over missing taxa in target tree
    for (int i = missingStart; i < missingEnd; i++) {
        int missingTaxon = missingTaxaFlat[i];
        
        // Check if this taxon is in the source range
        int posInSource = inverseIndex[sourceTree * numTaxa + missingTaxon];
        
        // If taxon exists in source tree AND is within range, count it
        if (posInSource >= 0 && posInSource >= rangeStart && posInSource < rangeEnd) {
            count++;
        }
    }
    
    return count;
}

// Legacy device function to calculate intersection cardinality using BitSets
__device__ int intersectionCardinality(unsigned long long* bits1, unsigned long long* bits2, int bitsetSize) {
    int count = 0;
    for (int i = 0; i < bitsetSize; i++) {
        unsigned long long intersection = bits1[i] & bits2[i];
        count += __popcll(intersection);  // Count set bits in 64-bit words
    }
    return count;
}

// Memory-optimized device function to calculate range intersection using inverse indices
__device__ int compactRangeIntersection(
    int tree1, int start1, int end1,
    int tree2, int start2, int end2,
    int* inverseIndex,    // Flattened [tree*numTaxa + taxon] = position
    int* orderings,       // Flattened [tree*numTaxa + position] = taxon
    int numTaxa
) {
    // Choose smaller range for iteration (complexity optimization)
    int size1 = end1 - start1;
    int size2 = end2 - start2;
    
    int smallTree, smallStart, smallEnd, largeTree, largeStart, largeEnd;
    if (size1 <= size2) {
        smallTree = tree1; smallStart = start1; smallEnd = end1;
        largeTree = tree2; largeStart = start2; largeEnd = end2;
    } else {
        smallTree = tree2; smallStart = start2; smallEnd = end2;
        largeTree = tree1; largeStart = start1; largeEnd = end1;
    }
    
    int count = 0;
    
    // Iterate over smaller range
    for (int pos = smallStart; pos < smallEnd; pos++) {
        // Get taxon ID at this position in small tree
        int taxonId = orderings[smallTree * numTaxa + pos];
        
        // Find position of this taxon in large tree using inverse index
        int positionInLargeTree = inverseIndex[largeTree * numTaxa + taxonId];
        
        // Check if taxon falls within large tree's range
        if (positionInLargeTree >= largeStart && positionInLargeTree < largeEnd) {
            count++;
        }
    }
    
    return count;
}

// Device function to calculate score for a pair of bipartitions
__device__ double calculateScore(Bipartition stb1, Bipartition stb2) {
    // First configuration: (A|B) with (X|Y)
    int p1 = intersectionCardinality(stb1.cluster1, stb2.cluster1, stb1.bitsetSize);
    int p2 = intersectionCardinality(stb1.cluster2, stb2.cluster2, stb1.bitsetSize);
    
    double score1 = 0.0;
    if (p1 + p2 >= 2) {
        score1 = p1 * p2 * (p1 + p2 - 2) / 2.0;
    }
    
    // Second configuration: (A|B) with (Y|X)
    p1 = intersectionCardinality(stb1.cluster1, stb2.cluster2, stb1.bitsetSize);
    p2 = intersectionCardinality(stb1.cluster2, stb2.cluster1, stb1.bitsetSize);
    
    double score2 = 0.0;
    if (p1 + p2 >= 2) {
        score2 = p1 * p2 * (p1 + p2 - 2) / 2.0;
    }
    
    return score1 + score2;
}

// Memory-optimized device function to calculate score using compact ranges
__device__ double calculateCompactScore(
    CompactBipartition bip1, CompactBipartition bip2,
    int* inverseIndex, int* orderings, int numTaxa
) {
    // Calculate four intersection sizes: AA, AB, BA, BB
    int aa = compactRangeIntersection(
        bip1.geneTreeIndex, bip1.leftStart, bip1.leftEnd,
        bip2.geneTreeIndex, bip2.leftStart, bip2.leftEnd,
        inverseIndex, orderings, numTaxa);
        
    int bb = compactRangeIntersection(
        bip1.geneTreeIndex, bip1.rightStart, bip1.rightEnd,
        bip2.geneTreeIndex, bip2.rightStart, bip2.rightEnd,
        inverseIndex, orderings, numTaxa);
        
    int ab = compactRangeIntersection(
        bip1.geneTreeIndex, bip1.leftStart, bip1.leftEnd,
        bip2.geneTreeIndex, bip2.rightStart, bip2.rightEnd,
        inverseIndex, orderings, numTaxa);
        
    int ba = compactRangeIntersection(
        bip1.geneTreeIndex, bip1.rightStart, bip1.rightEnd,
        bip2.geneTreeIndex, bip2.leftStart, bip2.leftEnd,
        inverseIndex, orderings, numTaxa);
    
    // Apply same scoring formula as original implementation
    double score1 = 0.0;
    if (aa + bb >= 2) {
        score1 = aa * bb * (aa + bb - 2) / 2.0;
    }
    
    double score2 = 0.0;
    if (ab + ba >= 2) {
        score2 = ab * ba * (ab + ba - 2) / 2.0;
    }
    
    return score1 + score2;
}

// Device function to calculate score for mixed bipartition (cross-tree recombination)
// Key difference: left and right sides use separate tree indices
__device__ double calculateMixedCompactScore(
    MixedCompactBipartition mixed, CompactBipartition geneTree,
    int* inverseIndex, int* orderings, int numTaxa
) {
    // Calculate four intersection sizes: AA, AB, BA, BB
    // Use leftTreeIndex for left side, rightTreeIndex for right side
    int aa = compactRangeIntersection(
        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
        geneTree.geneTreeIndex, geneTree.leftStart, geneTree.leftEnd,
        inverseIndex, orderings, numTaxa);
        
    int bb = compactRangeIntersection(
        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
        geneTree.geneTreeIndex, geneTree.rightStart, geneTree.rightEnd,
        inverseIndex, orderings, numTaxa);
        
    int ab = compactRangeIntersection(
        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
        geneTree.geneTreeIndex, geneTree.rightStart, geneTree.rightEnd,
        inverseIndex, orderings, numTaxa);
        
    int ba = compactRangeIntersection(
        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
        geneTree.geneTreeIndex, geneTree.leftStart, geneTree.leftEnd,
        inverseIndex, orderings, numTaxa);
    
    // Apply same scoring formula
    double score1 = 0.0;
    if (aa + bb >= 2) {
        score1 = aa * bb * (aa + bb - 2) / 2.0;
    }
    
    double score2 = 0.0;
    if (ab + ba >= 2) {
        score2 = ab * ba * (ab + ba - 2) / 2.0;
    }
    
    return score1 + score2;
}

// ============================================================================
// QUARTET SCORING FUNCTIONS (ASTRAL-style)
// ============================================================================

// Device function to calculate ASTRAL-style quartet score using compact ranges
// Uses tripartition (A|B|C) vs (X|Y|Z) with only 4 intersections + derivation
__device__ double calculateCompactQuartetScore(
    CompactBipartition candidate, CompactBipartition geneTree,
    int* inverseIndex, int* orderings,
    int* presentCount,      // presentCount[treeIndex] = |Sg|
    int* missingTaxaFlat,   // Flattened missing taxa list
    int* missingOffsets,    // Start offsets in missingTaxaFlat
    int* missingCounts,     // Number of missing taxa per tree
    int numTaxa
) {
    int geneTreeIdx = geneTree.geneTreeIndex;
    int sgSize = presentCount[geneTreeIdx];
    
    // Get sizes of gene tree bipartition sides
    int xSize = geneTree.leftEnd - geneTree.leftStart;
    int ySize = geneTree.rightEnd - geneTree.rightStart;
    int zSize = sgSize - xSize - ySize;
    if (zSize < 0) zSize = 0;
    
    // Get restricted sizes of candidate sides: a = |A ∩ Sg|, b = |B ∩ Sg|
    int aSize = (candidate.leftEnd - candidate.leftStart);
    int bSize = (candidate.rightEnd - candidate.rightStart);
    
    // Count missing taxa in each candidate side
    int aMissing = countMissingInRange(
        candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd,
        geneTreeIdx, inverseIndex, orderings,
        missingTaxaFlat, missingOffsets, missingCounts, numTaxa);
    int bMissing = countMissingInRange(
        candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd,
        geneTreeIdx, inverseIndex, orderings,
        missingTaxaFlat, missingOffsets, missingCounts, numTaxa);
    
    int a = aSize - aMissing;  // |A ∩ Sg|
    int b = bSize - bMissing;  // |B ∩ Sg|
    int c = sgSize - a - b;    // |C ∩ Sg| = |Sg| - a - b
    if (c < 0) c = 0;
    
    // Compute 4 intersections
    int ax = compactRangeIntersection(
        candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd,
        geneTreeIdx, geneTree.leftStart, geneTree.leftEnd,
        inverseIndex, orderings, numTaxa);
    
    int ay = compactRangeIntersection(
        candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd,
        geneTreeIdx, geneTree.rightStart, geneTree.rightEnd,
        inverseIndex, orderings, numTaxa);
    
    int bx = compactRangeIntersection(
        candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd,
        geneTreeIdx, geneTree.leftStart, geneTree.leftEnd,
        inverseIndex, orderings, numTaxa);
    
    int by = compactRangeIntersection(
        candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd,
        geneTreeIdx, geneTree.rightStart, geneTree.rightEnd,
        inverseIndex, orderings, numTaxa);
    
    // Derive remaining cells using row and column sums
    int az = a - ax - ay;
    int bz = b - bx - by;
    int cx = xSize - ax - bx;
    int cy = ySize - ay - by;
    int cz = zSize - az - bz;
    
    // Clamp negative values to 0 (can happen due to numerical issues)
    if (az < 0) az = 0;
    if (bz < 0) bz = 0;
    if (cx < 0) cx = 0;
    if (cy < 0) cy = 0;
    if (cz < 0) cz = 0;
    
    // Build flattened 3x3 grid and compute QI
    int grid[9] = {ax, ay, az, bx, by, bz, cx, cy, cz};
    
    return computeQuartetQI(grid);
}

// Device function to calculate quartet score for mixed bipartition
__device__ double calculateMixedCompactQuartetScore(
    MixedCompactBipartition mixed, CompactBipartition geneTree,
    int* inverseIndex, int* orderings,
    int* presentCount,
    int* missingTaxaFlat,
    int* missingOffsets,
    int* missingCounts,
    int numTaxa
) {
    int geneTreeIdx = geneTree.geneTreeIndex;
    int sgSize = presentCount[geneTreeIdx];
    
    // Get sizes of gene tree bipartition sides
    int xSize = geneTree.leftEnd - geneTree.leftStart;
    int ySize = geneTree.rightEnd - geneTree.rightStart;
    int zSize = sgSize - xSize - ySize;
    if (zSize < 0) zSize = 0;
    
    // Get restricted sizes - use respective tree indices for each side
    int aSize = mixed.leftEnd - mixed.leftStart;
    int bSize = mixed.rightEnd - mixed.rightStart;
    
    int aMissing = countMissingInRange(
        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
        geneTreeIdx, inverseIndex, orderings,
        missingTaxaFlat, missingOffsets, missingCounts, numTaxa);
    int bMissing = countMissingInRange(
        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
        geneTreeIdx, inverseIndex, orderings,
        missingTaxaFlat, missingOffsets, missingCounts, numTaxa);
    
    int a = aSize - aMissing;
    int b = bSize - bMissing;
    int c = sgSize - a - b;
    if (c < 0) c = 0;
    
    // Compute 4 intersections using respective tree indices
    int ax = compactRangeIntersection(
        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
        geneTreeIdx, geneTree.leftStart, geneTree.leftEnd,
        inverseIndex, orderings, numTaxa);
    
    int ay = compactRangeIntersection(
        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
        geneTreeIdx, geneTree.rightStart, geneTree.rightEnd,
        inverseIndex, orderings, numTaxa);
    
    int bx = compactRangeIntersection(
        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
        geneTreeIdx, geneTree.leftStart, geneTree.leftEnd,
        inverseIndex, orderings, numTaxa);
    
    int by = compactRangeIntersection(
        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
        geneTreeIdx, geneTree.rightStart, geneTree.rightEnd,
        inverseIndex, orderings, numTaxa);
    
    // Derive remaining cells
    int az = a - ax - ay;
    int bz = b - bx - by;
    int cx = xSize - ax - bx;
    int cy = ySize - ay - by;
    int cz = zSize - az - bz;
    
    // Clamp negative values
    if (az < 0) az = 0;
    if (bz < 0) bz = 0;
    if (cx < 0) cx = 0;
    if (cy < 0) cy = 0;
    if (cz < 0) cz = 0;
    
    int grid[9] = {ax, ay, az, bx, by, bz, cx, cy, cz};
    
    return computeQuartetQI(grid);
}

// Kernel to calculate weights for all candidate bipartitions
__global__ void calculateWeightsKernel(
    Bipartition* candidates,           // Array of candidate bipartitions
    Bipartition* geneTreeBips,         // Array of gene tree bipartitions
    int* frequencies,                  // Array of frequencies for gene tree bipartitions
    double* weights,                   // Output array for weights
    int numCandidates,                 // Number of candidate bipartitions
    int numGeneTreeBips                // Number of gene tree bipartitions
) {
    int candidateIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (candidateIdx >= numCandidates) return;
    
    double totalScore = 0.0;
    Bipartition candidate = candidates[candidateIdx];
    
    for (int i = 0; i < numGeneTreeBips; i++) {
        double score = calculateScore(candidate, geneTreeBips[i]);
        totalScore += score * frequencies[i];
    }
    
    weights[candidateIdx] = totalScore;
}

// Host function to launch the kernel
extern "C" {
    void launchWeightCalculation(
        Bipartition* hCandidates,
        Bipartition* hGeneTreeBips,
        int* hFrequencies,
        double* hWeights,
        int numCandidates,
        int numGeneTreeBips,
        int bitsetSize
    ) {
        // Device allocations
        Bipartition *dCandidates, *dGeneTreeBips;
        int *dFrequencies;
        double *dWeights;
        unsigned long long *dCandidateCluster1, *dCandidateCluster2;
        unsigned long long *dGeneTreeCluster1, *dGeneTreeCluster2;
        
        // Allocate device memory for Bipartition structures
        cudaMalloc(&dCandidates, numCandidates * sizeof(Bipartition));
        cudaMalloc(&dGeneTreeBips, numGeneTreeBips * sizeof(Bipartition));
        
        // Allocate device memory for bit arrays
        cudaMalloc(&dCandidateCluster1, numCandidates * bitsetSize * sizeof(unsigned long long));
        cudaMalloc(&dCandidateCluster2, numCandidates * bitsetSize * sizeof(unsigned long long));
        cudaMalloc(&dGeneTreeCluster1, numGeneTreeBips * bitsetSize * sizeof(unsigned long long));
        cudaMalloc(&dGeneTreeCluster2, numGeneTreeBips * bitsetSize * sizeof(unsigned long long));
        
        // Copy bit arrays to device
        for (int i = 0; i < numCandidates; i++) {
            cudaMemcpy(dCandidateCluster1 + i * bitsetSize, hCandidates[i].cluster1, bitsetSize * sizeof(unsigned long long), cudaMemcpyHostToDevice);
            cudaMemcpy(dCandidateCluster2 + i * bitsetSize, hCandidates[i].cluster2, bitsetSize * sizeof(unsigned long long), cudaMemcpyHostToDevice);
            
            Bipartition dBip;
            dBip.cluster1 = dCandidateCluster1 + i * bitsetSize;
            dBip.cluster2 = dCandidateCluster2 + i * bitsetSize;
            dBip.bitsetSize = bitsetSize;
            
            cudaMemcpy(&dCandidates[i], &dBip, sizeof(Bipartition), cudaMemcpyHostToDevice);
        }
        
        for (int i = 0; i < numGeneTreeBips; i++) {
            cudaMemcpy(dGeneTreeCluster1 + i * bitsetSize, hGeneTreeBips[i].cluster1, bitsetSize * sizeof(unsigned long long), cudaMemcpyHostToDevice);
            cudaMemcpy(dGeneTreeCluster2 + i * bitsetSize, hGeneTreeBips[i].cluster2, bitsetSize * sizeof(unsigned long long), cudaMemcpyHostToDevice);
            
            Bipartition dBip;
            dBip.cluster1 = dGeneTreeCluster1 + i * bitsetSize;
            dBip.cluster2 = dGeneTreeCluster2 + i * bitsetSize;
            dBip.bitsetSize = bitsetSize;
            
            cudaMemcpy(&dGeneTreeBips[i], &dBip, sizeof(Bipartition), cudaMemcpyHostToDevice);
        }
        
        // Allocate and copy other arrays
        cudaMalloc(&dFrequencies, numGeneTreeBips * sizeof(int));
        cudaMalloc(&dWeights, numCandidates * sizeof(double));
        cudaMemcpy(dFrequencies, hFrequencies, numGeneTreeBips * sizeof(int), cudaMemcpyHostToDevice);
        
        // Launch kernel
        int blockSize = 256;
        int gridSize = (numCandidates + blockSize - 1) / blockSize;
        
        calculateWeightsKernel<<<gridSize, blockSize>>>(
            dCandidates, dGeneTreeBips, dFrequencies, dWeights,
            numCandidates, numGeneTreeBips
        );
        
        // Copy results back to host
        cudaMemcpy(hWeights, dWeights, numCandidates * sizeof(double), cudaMemcpyDeviceToHost);
        
        // Cleanup
        cudaFree(dCandidates);
        cudaFree(dGeneTreeBips);
        cudaFree(dFrequencies);
        cudaFree(dWeights);
        cudaFree(dCandidateCluster1);
        cudaFree(dCandidateCluster2);
        cudaFree(dGeneTreeCluster1);
        cudaFree(dGeneTreeCluster2);
    }
    
    // Memory-optimized compact kernel for range-based bipartitions
    __global__ void calculateCompactWeightsKernel(
        CompactBipartition* candidates,    // Array of compact candidate bipartitions
        CompactBipartition* geneTreeBips,  // Array of compact gene tree bipartitions
        int* frequencies,                  // Array of frequencies for gene tree bipartitions
        double* weights,                   // Output array for weights
        int* inverseIndex,                 // Flattened inverse index [tree*numTaxa + taxon] = position
        int* orderings,                    // Flattened orderings [tree*numTaxa + position] = taxon
        int numCandidates,                 // Number of candidate bipartitions
        int numGeneTreeBips,               // Number of gene tree bipartitions
        int numTaxa                        // Number of taxa
    ) {
        int candidateIdx = blockIdx.x * blockDim.x + threadIdx.x;
        if (candidateIdx >= numCandidates) return;
        
        double totalScore = 0.0;
        CompactBipartition candidate = candidates[candidateIdx];
        
        // Calculate score against all gene tree bipartitions
        for (int i = 0; i < numGeneTreeBips; i++) {
            double score = calculateCompactScore(candidate, geneTreeBips[i], inverseIndex, orderings, numTaxa);
            totalScore += score * frequencies[i];
        }
        
        weights[candidateIdx] = totalScore;
    }
    
    // Host function to launch the compact kernel
    void launchCompactWeightCalculation(
        CompactBipartition* hCandidates,
        CompactBipartition* hGeneTreeBips,
        int* hFrequencies,
        double* hWeights,
        int* hInverseIndex,
        int* hOrderings,
        int numCandidates,
        int numGeneTreeBips,
        int numTrees,
        int numTaxa
    ) {
        printf("==== LAUNCHING COMPACT GPU KERNEL ====\n");
        printf("Candidates: %d, Gene tree bips: %d, Trees: %d, Taxa: %d\n", 
               numCandidates, numGeneTreeBips, numTrees, numTaxa);
        
        // Device allocations
        CompactBipartition *dCandidates, *dGeneTreeBips;
        int *dFrequencies, *dInverseIndex, *dOrderings;
        double *dWeights;
        
        // Calculate memory sizes
        size_t candidateSize = numCandidates * sizeof(CompactBipartition);
        size_t geneTreeSize = numGeneTreeBips * sizeof(CompactBipartition);
        size_t frequencySize = numGeneTreeBips * sizeof(int);
        size_t weightsSize = numCandidates * sizeof(double);
        size_t inverseIndexSize = (size_t)numTrees * numTaxa * sizeof(int);
        size_t orderingSize = (size_t)numTrees * numTaxa * sizeof(int);
        
        printf("Memory allocations:\n");
        printf("  Candidates: %zu KB\n", candidateSize / 1024);
        printf("  Gene trees: %zu KB\n", geneTreeSize / 1024);
        printf("  Inverse index: %zu MB\n", inverseIndexSize / (1024 * 1024));
        printf("  Orderings: %zu MB\n", orderingSize / (1024 * 1024));
        printf("  Total: %zu MB\n", (candidateSize + geneTreeSize + frequencySize + weightsSize + inverseIndexSize + orderingSize) / (1024 * 1024));
        
        // Allocate device memory
        cudaMalloc(&dCandidates, candidateSize);
        cudaMalloc(&dGeneTreeBips, geneTreeSize);
        cudaMalloc(&dFrequencies, frequencySize);
        cudaMalloc(&dWeights, weightsSize);
        cudaMalloc(&dInverseIndex, inverseIndexSize);
        cudaMalloc(&dOrderings, orderingSize);
        
        // Copy data to device
        printf("Copying data to GPU...\n");
        cudaMemcpy(dCandidates, hCandidates, candidateSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dGeneTreeBips, hGeneTreeBips, geneTreeSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dFrequencies, hFrequencies, frequencySize, cudaMemcpyHostToDevice);
        cudaMemcpy(dInverseIndex, hInverseIndex, inverseIndexSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dOrderings, hOrderings, orderingSize, cudaMemcpyHostToDevice);
        
        // Configure kernel launch parameters
        int blockSize = 256;  // Threads per block
        int gridSize = (numCandidates + blockSize - 1) / blockSize;  // Blocks per grid
        
        printf("Kernel configuration: %d blocks x %d threads = %d total threads\n", 
               gridSize, blockSize, gridSize * blockSize);
        
        // Launch kernel
        printf("Launching compact weight calculation kernel...\n");
        calculateCompactWeightsKernel<<<gridSize, blockSize>>>(
            dCandidates, dGeneTreeBips, dFrequencies, dWeights,
            dInverseIndex, dOrderings, numCandidates, numGeneTreeBips, numTaxa
        );
        
        // Check for kernel launch errors
        cudaError_t kernelError = cudaGetLastError();
        if (kernelError != cudaSuccess) {
            printf("CUDA kernel launch error: %s\n", cudaGetErrorString(kernelError));
            return;
        }
        
        // Wait for kernel to complete
        cudaDeviceSynchronize();
        
        // Check for kernel execution errors
        cudaError_t syncError = cudaGetLastError();
        if (syncError != cudaSuccess) {
            printf("CUDA kernel execution error: %s\n", cudaGetErrorString(syncError));
            return;
        }
        
        // Copy results back to host
        printf("Copying results back to host...\n");
        cudaMemcpy(hWeights, dWeights, weightsSize, cudaMemcpyDeviceToHost);
        
        // Free device memory
        cudaFree(dCandidates);
        cudaFree(dGeneTreeBips);
        cudaFree(dFrequencies);
        cudaFree(dWeights);
        cudaFree(dInverseIndex);
        cudaFree(dOrderings);
        
        printf("==== COMPACT GPU KERNEL COMPLETED SUCCESSFULLY ====\n");
    }
    
    // Mixed bipartition kernel for cross-tree recombination
    __global__ void calculateMixedCompactWeightsKernel(
        MixedCompactBipartition* mixedCandidates,  // Array of mixed bipartition candidates
        CompactBipartition* geneTreeBips,           // Array of compact gene tree bipartitions
        int* frequencies,                           // Array of frequencies for gene tree bipartitions
        double* weights,                            // Output array for weights
        int* inverseIndex,                          // Flattened inverse index
        int* orderings,                             // Flattened orderings
        int numCandidates,                          // Number of mixed bipartition candidates
        int numGeneTreeBips,                        // Number of gene tree bipartitions
        int numTaxa                                 // Number of taxa
    ) {
        int candidateIdx = blockIdx.x * blockDim.x + threadIdx.x;
        if (candidateIdx >= numCandidates) return;
        
        double totalScore = 0.0;
        MixedCompactBipartition candidate = mixedCandidates[candidateIdx];
        
        // Calculate score against all gene tree bipartitions
        for (int i = 0; i < numGeneTreeBips; i++) {
            double score = calculateMixedCompactScore(candidate, geneTreeBips[i], inverseIndex, orderings, numTaxa);
            totalScore += score * frequencies[i];
        }
        
        weights[candidateIdx] = totalScore;
    }
    
    // Host function to launch the mixed bipartition kernel
    void launchMixedWeightCalculation(
        MixedCompactBipartition* hMixedCandidates,
        CompactBipartition* hGeneTreeBips,
        int* hFrequencies,
        double* hWeights,
        int* hInverseIndex,
        int* hOrderings,
        int numCandidates,
        int numGeneTreeBips,
        int numTrees,
        int numTaxa
    ) {
        printf("==== LAUNCHING MIXED BIPARTITION GPU KERNEL ====\n");
        printf("Mixed candidates: %d, Gene tree bips: %d, Trees: %d, Taxa: %d\n", 
               numCandidates, numGeneTreeBips, numTrees, numTaxa);
        
        // Device allocations
        MixedCompactBipartition *dMixedCandidates;
        CompactBipartition *dGeneTreeBips;
        int *dFrequencies, *dInverseIndex, *dOrderings;
        double *dWeights;
        
        // Calculate memory sizes
        size_t mixedCandidateSize = numCandidates * sizeof(MixedCompactBipartition);
        size_t geneTreeSize = numGeneTreeBips * sizeof(CompactBipartition);
        size_t frequencySize = numGeneTreeBips * sizeof(int);
        size_t weightsSize = numCandidates * sizeof(double);
        size_t inverseIndexSize = (size_t)numTrees * numTaxa * sizeof(int);
        size_t orderingSize = (size_t)numTrees * numTaxa * sizeof(int);
        
        printf("Memory allocations:\n");
        printf("  Mixed candidates: %zu KB\n", mixedCandidateSize / 1024);
        printf("  Gene trees: %zu KB\n", geneTreeSize / 1024);
        printf("  Inverse index: %zu MB\n", inverseIndexSize / (1024 * 1024));
        printf("  Orderings: %zu MB\n", orderingSize / (1024 * 1024));
        printf("  Total: %zu MB\n", (mixedCandidateSize + geneTreeSize + frequencySize + weightsSize + inverseIndexSize + orderingSize) / (1024 * 1024));
        
        // Allocate device memory
        cudaMalloc(&dMixedCandidates, mixedCandidateSize);
        cudaMalloc(&dGeneTreeBips, geneTreeSize);
        cudaMalloc(&dFrequencies, frequencySize);
        cudaMalloc(&dWeights, weightsSize);
        cudaMalloc(&dInverseIndex, inverseIndexSize);
        cudaMalloc(&dOrderings, orderingSize);
        
        // Copy data to device
        printf("Copying data to GPU...\n");
        cudaMemcpy(dMixedCandidates, hMixedCandidates, mixedCandidateSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dGeneTreeBips, hGeneTreeBips, geneTreeSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dFrequencies, hFrequencies, frequencySize, cudaMemcpyHostToDevice);
        cudaMemcpy(dInverseIndex, hInverseIndex, inverseIndexSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dOrderings, hOrderings, orderingSize, cudaMemcpyHostToDevice);
        
        // Configure kernel launch parameters
        int blockSize = 256;
        int gridSize = (numCandidates + blockSize - 1) / blockSize;
        
        printf("Kernel configuration: %d blocks x %d threads = %d total threads\n", 
               gridSize, blockSize, gridSize * blockSize);
        
        // Launch kernel
        printf("Launching mixed bipartition weight calculation kernel...\n");
        calculateMixedCompactWeightsKernel<<<gridSize, blockSize>>>(
            dMixedCandidates, dGeneTreeBips, dFrequencies, dWeights,
            dInverseIndex, dOrderings, numCandidates, numGeneTreeBips, numTaxa
        );
        
        // Check for kernel launch errors
        cudaError_t kernelError = cudaGetLastError();
        if (kernelError != cudaSuccess) {
            printf("CUDA kernel launch error: %s\n", cudaGetErrorString(kernelError));
            return;
        }
        
        // Wait for kernel to complete
        cudaDeviceSynchronize();
        
        // Check for kernel execution errors
        cudaError_t syncError = cudaGetLastError();
        if (syncError != cudaSuccess) {
            printf("CUDA kernel execution error: %s\n", cudaGetErrorString(syncError));
            return;
        }
        
        // Copy results back to host
        printf("Copying results back to host...\n");
        cudaMemcpy(hWeights, dWeights, weightsSize, cudaMemcpyDeviceToHost);
        
        // Free device memory
        cudaFree(dMixedCandidates);
        cudaFree(dGeneTreeBips);
        cudaFree(dFrequencies);
        cudaFree(dWeights);
        cudaFree(dInverseIndex);
        cudaFree(dOrderings);
        
        printf("==== MIXED BIPARTITION GPU KERNEL COMPLETED SUCCESSFULLY ====\n");
    }
    
    // ============================================================================
    // QUARTET SCORING GPU KERNELS (ASTRAL-style)
    // ============================================================================
    
    // Kernel for quartet weight calculation using compact bipartitions
    __global__ void calculateCompactQuartetWeightsKernel(
        CompactBipartition* candidates,
        CompactBipartition* geneTreeBips,
        int* frequencies,
        double* weights,
        int* inverseIndex,
        int* orderings,
        int* presentCount,        // presentCount[treeIndex] = |Sg|
        int* missingTaxaFlat,     // Flattened missing taxa list
        int* missingOffsets,      // Start offsets in missingTaxaFlat
        int* missingCounts,       // Number of missing taxa per tree
        int numCandidates,
        int numGeneTreeBips,
        int numTaxa
    ) {
        int candidateIdx = blockIdx.x * blockDim.x + threadIdx.x;
        if (candidateIdx >= numCandidates) return;
        
        double totalScore = 0.0;
        CompactBipartition candidate = candidates[candidateIdx];
        
        for (int i = 0; i < numGeneTreeBips; i++) {
            double score = calculateCompactQuartetScore(
                candidate, geneTreeBips[i],
                inverseIndex, orderings,
                presentCount, missingTaxaFlat, missingOffsets, missingCounts,
                numTaxa);
            totalScore += score * frequencies[i];
        }
        
        weights[candidateIdx] = totalScore;
    }
    
    // Kernel for quartet weight calculation with mixed bipartitions
    __global__ void calculateMixedCompactQuartetWeightsKernel(
        MixedCompactBipartition* mixedCandidates,
        CompactBipartition* geneTreeBips,
        int* frequencies,
        double* weights,
        int* inverseIndex,
        int* orderings,
        int* presentCount,
        int* missingTaxaFlat,
        int* missingOffsets,
        int* missingCounts,
        int numCandidates,
        int numGeneTreeBips,
        int numTaxa
    ) {
        int candidateIdx = blockIdx.x * blockDim.x + threadIdx.x;
        if (candidateIdx >= numCandidates) return;
        
        double totalScore = 0.0;
        MixedCompactBipartition candidate = mixedCandidates[candidateIdx];
        
        for (int i = 0; i < numGeneTreeBips; i++) {
            double score = calculateMixedCompactQuartetScore(
                candidate, geneTreeBips[i],
                inverseIndex, orderings,
                presentCount, missingTaxaFlat, missingOffsets, missingCounts,
                numTaxa);
            totalScore += score * frequencies[i];
        }
        
        weights[candidateIdx] = totalScore;
    }
    
    // Host function to launch quartet weight calculation
    void launchQuartetWeightCalculation(
        CompactBipartition* hCandidates,
        CompactBipartition* hGeneTreeBips,
        int* hFrequencies,
        double* hWeights,
        int* hInverseIndex,
        int* hOrderings,
        int* hPresentCount,       // presentCount[treeIndex] = |Sg|
        int* hMissingTaxaFlat,    // Flattened missing taxa list
        int* hMissingOffsets,     // Start offsets in missingTaxaFlat
        int* hMissingCounts,      // Number of missing taxa per tree
        int numCandidates,
        int numGeneTreeBips,
        int numTrees,
        int numTaxa,
        int totalMissingTaxa      // Total number of entries in missingTaxaFlat
    ) {
        printf("==== LAUNCHING QUARTET SCORING GPU KERNEL ====\n");
        printf("Candidates: %d, Gene tree bips: %d, Trees: %d, Taxa: %d\n",
               numCandidates, numGeneTreeBips, numTrees, numTaxa);
        printf("Total missing taxa entries: %d\n", totalMissingTaxa);
        
        // Device allocations
        CompactBipartition *dCandidates, *dGeneTreeBips;
        int *dFrequencies, *dInverseIndex, *dOrderings;
        int *dPresentCount, *dMissingTaxaFlat, *dMissingOffsets, *dMissingCounts;
        double *dWeights;
        
        // Calculate memory sizes
        size_t candidateSize = numCandidates * sizeof(CompactBipartition);
        size_t geneTreeSize = numGeneTreeBips * sizeof(CompactBipartition);
        size_t frequencySize = numGeneTreeBips * sizeof(int);
        size_t weightsSize = numCandidates * sizeof(double);
        size_t inverseIndexSize = (size_t)numTrees * numTaxa * sizeof(int);
        size_t orderingSize = (size_t)numTrees * numTaxa * sizeof(int);
        size_t presentCountSize = numTrees * sizeof(int);
        size_t missingTaxaFlatSize = totalMissingTaxa * sizeof(int);
        size_t missingOffsetsSize = numTrees * sizeof(int);
        size_t missingCountsSize = numTrees * sizeof(int);
        
        printf("Memory allocations:\n");
        printf("  Candidates: %zu KB\n", candidateSize / 1024);
        printf("  Gene trees: %zu KB\n", geneTreeSize / 1024);
        printf("  Inverse index: %zu MB\n", inverseIndexSize / (1024 * 1024));
        printf("  Missing taxa data: %zu KB\n", (presentCountSize + missingTaxaFlatSize + missingOffsetsSize + missingCountsSize) / 1024);
        
        // Allocate device memory
        cudaMalloc(&dCandidates, candidateSize);
        cudaMalloc(&dGeneTreeBips, geneTreeSize);
        cudaMalloc(&dFrequencies, frequencySize);
        cudaMalloc(&dWeights, weightsSize);
        cudaMalloc(&dInverseIndex, inverseIndexSize);
        cudaMalloc(&dOrderings, orderingSize);
        cudaMalloc(&dPresentCount, presentCountSize);
        cudaMalloc(&dMissingTaxaFlat, missingTaxaFlatSize > 0 ? missingTaxaFlatSize : sizeof(int));
        cudaMalloc(&dMissingOffsets, missingOffsetsSize);
        cudaMalloc(&dMissingCounts, missingCountsSize);
        
        // Copy data to device
        printf("Copying data to GPU...\n");
        cudaMemcpy(dCandidates, hCandidates, candidateSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dGeneTreeBips, hGeneTreeBips, geneTreeSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dFrequencies, hFrequencies, frequencySize, cudaMemcpyHostToDevice);
        cudaMemcpy(dInverseIndex, hInverseIndex, inverseIndexSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dOrderings, hOrderings, orderingSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dPresentCount, hPresentCount, presentCountSize, cudaMemcpyHostToDevice);
        if (totalMissingTaxa > 0) {
            cudaMemcpy(dMissingTaxaFlat, hMissingTaxaFlat, missingTaxaFlatSize, cudaMemcpyHostToDevice);
        }
        cudaMemcpy(dMissingOffsets, hMissingOffsets, missingOffsetsSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dMissingCounts, hMissingCounts, missingCountsSize, cudaMemcpyHostToDevice);
        
        // Configure kernel launch parameters
        int blockSize = 256;
        int gridSize = (numCandidates + blockSize - 1) / blockSize;
        
        printf("Kernel configuration: %d blocks x %d threads\n", gridSize, blockSize);
        
        // Launch kernel
        printf("Launching quartet weight calculation kernel...\n");
        calculateCompactQuartetWeightsKernel<<<gridSize, blockSize>>>(
            dCandidates, dGeneTreeBips, dFrequencies, dWeights,
            dInverseIndex, dOrderings,
            dPresentCount, dMissingTaxaFlat, dMissingOffsets, dMissingCounts,
            numCandidates, numGeneTreeBips, numTaxa
        );
        
        // Check for errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA kernel launch error: %s\n", cudaGetErrorString(err));
            return;
        }
        
        cudaDeviceSynchronize();
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA kernel execution error: %s\n", cudaGetErrorString(err));
            return;
        }
        
        // Copy results back
        cudaMemcpy(hWeights, dWeights, weightsSize, cudaMemcpyDeviceToHost);
        
        // Free device memory
        cudaFree(dCandidates);
        cudaFree(dGeneTreeBips);
        cudaFree(dFrequencies);
        cudaFree(dWeights);
        cudaFree(dInverseIndex);
        cudaFree(dOrderings);
        cudaFree(dPresentCount);
        cudaFree(dMissingTaxaFlat);
        cudaFree(dMissingOffsets);
        cudaFree(dMissingCounts);
        
        printf("==== QUARTET SCORING GPU KERNEL COMPLETED SUCCESSFULLY ====\n");
    }
    
    // Host function to launch mixed quartet weight calculation
    void launchMixedQuartetWeightCalculation(
        MixedCompactBipartition* hMixedCandidates,
        CompactBipartition* hGeneTreeBips,
        int* hFrequencies,
        double* hWeights,
        int* hInverseIndex,
        int* hOrderings,
        int* hPresentCount,
        int* hMissingTaxaFlat,
        int* hMissingOffsets,
        int* hMissingCounts,
        int numCandidates,
        int numGeneTreeBips,
        int numTrees,
        int numTaxa,
        int totalMissingTaxa
    ) {
        printf("==== LAUNCHING MIXED QUARTET SCORING GPU KERNEL ====\n");
        printf("Mixed candidates: %d, Gene tree bips: %d, Trees: %d, Taxa: %d\n",
               numCandidates, numGeneTreeBips, numTrees, numTaxa);
        
        // Device allocations
        MixedCompactBipartition *dMixedCandidates;
        CompactBipartition *dGeneTreeBips;
        int *dFrequencies, *dInverseIndex, *dOrderings;
        int *dPresentCount, *dMissingTaxaFlat, *dMissingOffsets, *dMissingCounts;
        double *dWeights;
        
        // Calculate memory sizes
        size_t mixedCandidateSize = numCandidates * sizeof(MixedCompactBipartition);
        size_t geneTreeSize = numGeneTreeBips * sizeof(CompactBipartition);
        size_t frequencySize = numGeneTreeBips * sizeof(int);
        size_t weightsSize = numCandidates * sizeof(double);
        size_t inverseIndexSize = (size_t)numTrees * numTaxa * sizeof(int);
        size_t orderingSize = (size_t)numTrees * numTaxa * sizeof(int);
        size_t presentCountSize = numTrees * sizeof(int);
        size_t missingTaxaFlatSize = totalMissingTaxa * sizeof(int);
        size_t missingOffsetsSize = numTrees * sizeof(int);
        size_t missingCountsSize = numTrees * sizeof(int);
        
        // Allocate device memory
        cudaMalloc(&dMixedCandidates, mixedCandidateSize);
        cudaMalloc(&dGeneTreeBips, geneTreeSize);
        cudaMalloc(&dFrequencies, frequencySize);
        cudaMalloc(&dWeights, weightsSize);
        cudaMalloc(&dInverseIndex, inverseIndexSize);
        cudaMalloc(&dOrderings, orderingSize);
        cudaMalloc(&dPresentCount, presentCountSize);
        cudaMalloc(&dMissingTaxaFlat, missingTaxaFlatSize > 0 ? missingTaxaFlatSize : sizeof(int));
        cudaMalloc(&dMissingOffsets, missingOffsetsSize);
        cudaMalloc(&dMissingCounts, missingCountsSize);
        
        // Copy data to device
        printf("Copying data to GPU...\n");
        cudaMemcpy(dMixedCandidates, hMixedCandidates, mixedCandidateSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dGeneTreeBips, hGeneTreeBips, geneTreeSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dFrequencies, hFrequencies, frequencySize, cudaMemcpyHostToDevice);
        cudaMemcpy(dInverseIndex, hInverseIndex, inverseIndexSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dOrderings, hOrderings, orderingSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dPresentCount, hPresentCount, presentCountSize, cudaMemcpyHostToDevice);
        if (totalMissingTaxa > 0) {
            cudaMemcpy(dMissingTaxaFlat, hMissingTaxaFlat, missingTaxaFlatSize, cudaMemcpyHostToDevice);
        }
        cudaMemcpy(dMissingOffsets, hMissingOffsets, missingOffsetsSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dMissingCounts, hMissingCounts, missingCountsSize, cudaMemcpyHostToDevice);
        
        // Configure and launch kernel
        int blockSize = 256;
        int gridSize = (numCandidates + blockSize - 1) / blockSize;
        
        printf("Launching mixed quartet weight calculation kernel...\n");
        calculateMixedCompactQuartetWeightsKernel<<<gridSize, blockSize>>>(
            dMixedCandidates, dGeneTreeBips, dFrequencies, dWeights,
            dInverseIndex, dOrderings,
            dPresentCount, dMissingTaxaFlat, dMissingOffsets, dMissingCounts,
            numCandidates, numGeneTreeBips, numTaxa
        );
        
        // Check for errors and copy results
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
        }
        cudaDeviceSynchronize();
        
        cudaMemcpy(hWeights, dWeights, weightsSize, cudaMemcpyDeviceToHost);
        
        // Free device memory
        cudaFree(dMixedCandidates);
        cudaFree(dGeneTreeBips);
        cudaFree(dFrequencies);
        cudaFree(dWeights);
        cudaFree(dInverseIndex);
        cudaFree(dOrderings);
        cudaFree(dPresentCount);
        cudaFree(dMissingTaxaFlat);
        cudaFree(dMissingOffsets);
        cudaFree(dMissingCounts);
        
        printf("==== MIXED QUARTET SCORING GPU KERNEL COMPLETED ====\n");
    }
} 