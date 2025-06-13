#include <stdio.h>
#include <cuda_runtime.h>

// Structure to represent a bipartition
struct Bipartition {
    unsigned long long* cluster1;  // First cluster bits
    unsigned long long* cluster2;  // Second cluster bits
    int bitsetSize;               // Size of bitset in 64-bit words
};

// Device function to calculate intersection cardinality
__device__ int intersectionCardinality(unsigned long long* bits1, unsigned long long* bits2, int bitsetSize) {
    int count = 0;
    for (int i = 0; i < bitsetSize; i++) {
        unsigned long long intersection = bits1[i] & bits2[i];
        count += __popcll(intersection);  // Count set bits in 64-bit word
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
} 