#include <stdio.h>
#include <cuda_runtime.h>

// Structure to represent a compact bipartition using ranges
struct CompactBipartition {
    int geneTreeIndex;    // Which gene tree this bipartition belongs to
    int leftStart;        // Start index of left cluster range (inclusive)
    int leftEnd;          // End index of left cluster range (exclusive)
    int rightStart;       // Start index of right cluster range (inclusive)
    int rightEnd;         // End index of right cluster range (exclusive)
};

// Device function to calculate range intersection using inverse index
__device__ int calculateRangeIntersection(
    int rangeA_start, int rangeA_end, int geneTreeA,
    int rangeB_start, int rangeB_end, int geneTreeB,
    int* inverseIndex,    // Flattened array: [geneTree * maxTaxa + taxonId] = position
    int* geneTreeOrderings, // Flattened array: [geneTree * maxTaxa + position] = taxonId
    int maxTaxa
) {
    int sizeA = rangeA_end - rangeA_start;
    int sizeB = rangeB_end - rangeB_start;
    
    // Iterate over the smaller range for efficiency
    if (sizeA <= sizeB) {
        int count = 0;
        for (int pos = rangeA_start; pos < rangeA_end; pos++) {
            int taxonId = geneTreeOrderings[geneTreeA * maxTaxa + pos];
            int positionInTarget = inverseIndex[geneTreeB * maxTaxa + taxonId];
            
            if (positionInTarget >= rangeB_start && positionInTarget < rangeB_end) {
                count++;
            }
        }
        return count;
    } else {
        int count = 0;
        for (int pos = rangeB_start; pos < rangeB_end; pos++) {
            int taxonId = geneTreeOrderings[geneTreeB * maxTaxa + pos];
            int positionInTarget = inverseIndex[geneTreeA * maxTaxa + taxonId];
            
            if (positionInTarget >= rangeA_start && positionInTarget < rangeA_end) {
                count++;
            }
        }
        return count;
    }
}

// Device function to calculate all four intersection counts for two compact bipartitions
__device__ void calculateAllIntersections(
    CompactBipartition bip1, CompactBipartition bip2,
    int* intersections, // Output array of size 4
    int* inverseIndex, int* geneTreeOrderings, int maxTaxa
) {
    // |A1 ∩ A2| - left cluster of bip1 with left cluster of bip2
    intersections[0] = calculateRangeIntersection(
        bip1.leftStart, bip1.leftEnd, bip1.geneTreeIndex,
        bip2.leftStart, bip2.leftEnd, bip2.geneTreeIndex,
        inverseIndex, geneTreeOrderings, maxTaxa
    );
    
    // |A1 ∩ B2| - left cluster of bip1 with right cluster of bip2
    intersections[1] = calculateRangeIntersection(
        bip1.leftStart, bip1.leftEnd, bip1.geneTreeIndex,
        bip2.rightStart, bip2.rightEnd, bip2.geneTreeIndex,
        inverseIndex, geneTreeOrderings, maxTaxa
    );
    
    // |B1 ∩ A2| - right cluster of bip1 with left cluster of bip2
    intersections[2] = calculateRangeIntersection(
        bip1.rightStart, bip1.rightEnd, bip1.geneTreeIndex,
        bip2.leftStart, bip2.leftEnd, bip2.geneTreeIndex,
        inverseIndex, geneTreeOrderings, maxTaxa
    );
    
    // |B1 ∩ B2| - right cluster of bip1 with right cluster of bip2
    intersections[3] = calculateRangeIntersection(
        bip1.rightStart, bip1.rightEnd, bip1.geneTreeIndex,
        bip2.rightStart, bip2.rightEnd, bip2.geneTreeIndex,
        inverseIndex, geneTreeOrderings, maxTaxa
    );
}

// Device function to calculate score for a pair of compact bipartitions
__device__ double calculateCompactScore(
    CompactBipartition bip1, CompactBipartition bip2,
    int* inverseIndex, int* geneTreeOrderings, int maxTaxa
) {
    int intersections[4];
    calculateAllIntersections(bip1, bip2, intersections, inverseIndex, geneTreeOrderings, maxTaxa);
    
    // First configuration: (A1|B1) with (A2|B2)
    int p1 = intersections[0]; // |A1 ∩ A2|
    int p2 = intersections[3]; // |B1 ∩ B2|
    
    double score1 = 0.0;
    if (p1 + p2 >= 2) {
        score1 = p1 * p2 * (p1 + p2 - 2) / 2.0;
    }
    
    // Second configuration: (A1|B1) with (B2|A2) - cross configuration
    p1 = intersections[1]; // |A1 ∩ B2|
    p2 = intersections[2]; // |B1 ∩ A2|
    
    double score2 = 0.0;
    if (p1 + p2 >= 2) {
        score2 = p1 * p2 * (p1 + p2 - 2) / 2.0;
    }
    
    return score1 + score2;
}

// Kernel to calculate weights for all candidate compact bipartitions
__global__ void calculateCompactWeightsKernel(
    CompactBipartition* candidates,        // Array of candidate compact bipartitions
    CompactBipartition* geneTreeBips,      // Array of gene tree compact bipartitions
    int* frequencies,                      // Array of frequencies for gene tree bipartitions
    double* weights,                       // Output array for weights
    int* inverseIndex,                     // Inverse index matrix (flattened)
    int* geneTreeOrderings,                // Gene tree orderings (flattened)
    int numCandidates,                     // Number of candidate bipartitions
    int numGeneTreeBips,                   // Number of gene tree bipartitions
    int maxTaxa                            // Maximum number of taxa
) {
    int candidateIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (candidateIdx >= numCandidates) return;
    
    double totalScore = 0.0;
    CompactBipartition candidate = candidates[candidateIdx];
    
    for (int i = 0; i < numGeneTreeBips; i++) {
        double score = calculateCompactScore(candidate, geneTreeBips[i], 
                                           inverseIndex, geneTreeOrderings, maxTaxa);
        totalScore += score * frequencies[i];
    }
    
    weights[candidateIdx] = totalScore;
}

// Host function to launch the memory-efficient kernel
extern "C" {
    void launchMemoryEfficientWeightCalculation(
        CompactBipartition* hCandidates,
        CompactBipartition* hGeneTreeBips,
        int* hFrequencies,
        double* hWeights,
        int* hInverseIndex,
        int* hGeneTreeOrderings,
        int numCandidates,
        int numGeneTreeBips,
        int maxTaxa,
        int numGeneTrees
    ) {
        // Device allocations
        CompactBipartition *dCandidates, *dGeneTreeBips;
        int *dFrequencies, *dInverseIndex, *dGeneTreeOrderings;
        double *dWeights;
        
        // Calculate memory sizes
        size_t candidatesSize = numCandidates * sizeof(CompactBipartition);
        size_t geneTreeBipsSize = numGeneTreeBips * sizeof(CompactBipartition);
        size_t frequenciesSize = numGeneTreeBips * sizeof(int);
        size_t weightsSize = numCandidates * sizeof(double);
        size_t inverseIndexSize = numGeneTrees * maxTaxa * sizeof(int);
        size_t orderingsSize = numGeneTrees * maxTaxa * sizeof(int);
        
        // Allocate device memory
        cudaMalloc(&dCandidates, candidatesSize);
        cudaMalloc(&dGeneTreeBips, geneTreeBipsSize);
        cudaMalloc(&dFrequencies, frequenciesSize);
        cudaMalloc(&dWeights, weightsSize);
        cudaMalloc(&dInverseIndex, inverseIndexSize);
        cudaMalloc(&dGeneTreeOrderings, orderingsSize);
        
        // Copy data to device
        cudaMemcpy(dCandidates, hCandidates, candidatesSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dGeneTreeBips, hGeneTreeBips, geneTreeBipsSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dFrequencies, hFrequencies, frequenciesSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dInverseIndex, hInverseIndex, inverseIndexSize, cudaMemcpyHostToDevice);
        cudaMemcpy(dGeneTreeOrderings, hGeneTreeOrderings, orderingsSize, cudaMemcpyHostToDevice);
        
        // Launch kernel
        int blockSize = 256;
        int gridSize = (numCandidates + blockSize - 1) / blockSize;
        
        calculateCompactWeightsKernel<<<gridSize, blockSize>>>(
            dCandidates, dGeneTreeBips, dFrequencies, dWeights,
            dInverseIndex, dGeneTreeOrderings,
            numCandidates, numGeneTreeBips, maxTaxa
        );
        
        // Copy results back to host
        cudaMemcpy(hWeights, dWeights, weightsSize, cudaMemcpyDeviceToHost);
        
        // Cleanup
        cudaFree(dCandidates);
        cudaFree(dGeneTreeBips);
        cudaFree(dFrequencies);
        cudaFree(dWeights);
        cudaFree(dInverseIndex);
        cudaFree(dGeneTreeOrderings);
    }
}
