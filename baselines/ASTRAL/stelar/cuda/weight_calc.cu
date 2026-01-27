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
} 