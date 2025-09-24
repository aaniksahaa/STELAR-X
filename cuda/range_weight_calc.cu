#include <stdio.h>
#include <cuda_runtime.h>
#include <algorithm>

// Structure to represent a range-based bipartition
struct RangeBipartition {
    int geneTreeIndex;  // Index of the gene tree
    int startPos;       // Start position in taxa ordering
    int endPos;         // End position (exclusive) in taxa ordering
};

// Device function to get taxa array for a range bipartition
__device__ void getRangeTaxa(const RangeBipartition& bip, const int* geneTreeOrderings, 
                            const int* orderingOffsets, const int* orderingLengths, 
                            int* result, int* resultSize) {
    if (bip.geneTreeIndex < 0 || bip.startPos < 0 || bip.endPos <= bip.startPos) {
        *resultSize = 0;
        return;
    }
    
    int offset = orderingOffsets[bip.geneTreeIndex];
    int length = orderingLengths[bip.geneTreeIndex];
    
    if (bip.endPos > length) {
        *resultSize = 0;
        return;
    }
    
    *resultSize = bip.endPos - bip.startPos;
    for (int i = 0; i < *resultSize; i++) {
        result[i] = geneTreeOrderings[offset + bip.startPos + i];
    }
}

// Device function to get complement taxa array for a range bipartition
__device__ void getRangeComplementTaxa(const RangeBipartition& bip, const int* geneTreeOrderings,
                                      const int* orderingOffsets, const int* orderingLengths,
                                      int* result, int* resultSize) {
    if (bip.geneTreeIndex < 0) {
        *resultSize = 0;
        return;
    }
    
    int offset = orderingOffsets[bip.geneTreeIndex];
    int length = orderingLengths[bip.geneTreeIndex];
    
    int idx = 0;
    // Add taxa before startPos
    for (int i = 0; i < bip.startPos && i < length; i++) {
        result[idx++] = geneTreeOrderings[offset + i];
    }
    // Add taxa after endPos
    for (int i = bip.endPos; i < length; i++) {
        result[idx++] = geneTreeOrderings[offset + i];
    }
    
    *resultSize = idx;
}

// Device function to calculate intersection size between two sorted arrays
__device__ int calculateIntersectionSize(int* arr1, int size1, int* arr2, int size2) {
    // Sort both arrays (bubble sort for small arrays in GPU)
    for (int i = 0; i < size1 - 1; i++) {
        for (int j = 0; j < size1 - i - 1; j++) {
            if (arr1[j] > arr1[j + 1]) {
                int temp = arr1[j];
                arr1[j] = arr1[j + 1];
                arr1[j + 1] = temp;
            }
        }
    }
    
    for (int i = 0; i < size2 - 1; i++) {
        for (int j = 0; j < size2 - i - 1; j++) {
            if (arr2[j] > arr2[j + 1]) {
                int temp = arr2[j];
                arr2[j] = arr2[j + 1];
                arr2[j + 1] = temp;
            }
        }
    }
    
    // Count intersection using two pointers
    int count = 0;
    int i = 0, j = 0;
    
    while (i < size1 && j < size2) {
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

// Device function to calculate score for a pair of range bipartitions
__device__ double calculateRangeScore(const RangeBipartition& stb1, const RangeBipartition& stb2,
                                     const int* geneTreeOrderings, const int* orderingOffsets,
                                     const int* orderingLengths, int maxTaxaCount) {
    
    // Allocate temporary arrays (shared memory would be better for large datasets)
    int* stb1_cluster1 = new int[maxTaxaCount];
    int* stb1_cluster2 = new int[maxTaxaCount];
    int* stb2_cluster1 = new int[maxTaxaCount];
    int* stb2_cluster2 = new int[maxTaxaCount];
    
    int stb1_cluster1_size, stb1_cluster2_size, stb2_cluster1_size, stb2_cluster2_size;
    
    // Get taxa arrays for both bipartitions
    getRangeTaxa(stb1, geneTreeOrderings, orderingOffsets, orderingLengths, stb1_cluster1, &stb1_cluster1_size);
    getRangeComplementTaxa(stb1, geneTreeOrderings, orderingOffsets, orderingLengths, stb1_cluster2, &stb1_cluster2_size);
    getRangeTaxa(stb2, geneTreeOrderings, orderingOffsets, orderingLengths, stb2_cluster1, &stb2_cluster1_size);
    getRangeComplementTaxa(stb2, geneTreeOrderings, orderingOffsets, orderingLengths, stb2_cluster2, &stb2_cluster2_size);
    
    // First configuration: (A|B) with (X|Y)
    int p1 = calculateIntersectionSize(stb1_cluster1, stb1_cluster1_size, stb2_cluster1, stb2_cluster1_size);
    int p2 = calculateIntersectionSize(stb1_cluster2, stb1_cluster2_size, stb2_cluster2, stb2_cluster2_size);
    
    double score1 = 0.0;
    if (p1 + p2 >= 2) {
        score1 = p1 * p2 * (p1 + p2 - 2) / 2.0;
    }
    
    // Second configuration: (A|B) with (Y|X) - cross configuration
    p1 = calculateIntersectionSize(stb1_cluster1, stb1_cluster1_size, stb2_cluster2, stb2_cluster2_size);
    p2 = calculateIntersectionSize(stb1_cluster2, stb1_cluster2_size, stb2_cluster1, stb2_cluster1_size);
    
    double score2 = 0.0;
    if (p1 + p2 >= 2) {
        score2 = p1 * p2 * (p1 + p2 - 2) / 2.0;
    }
    
    // Cleanup
    delete[] stb1_cluster1;
    delete[] stb1_cluster2;
    delete[] stb2_cluster1;
    delete[] stb2_cluster2;
    
    return score1 + score2;
}

// Kernel to calculate weights for all candidate range bipartitions
__global__ void calculateRangeWeightsKernel(
    const RangeBipartition* candidates,      // Array of candidate bipartitions
    const RangeBipartition* geneTreeBips,    // Array of gene tree bipartitions
    const int* frequencies,                  // Array of frequencies for gene tree bipartitions
    double* weights,                         // Output array for weights
    const int* geneTreeOrderings,            // Flattened array of all gene tree orderings
    const int* orderingOffsets,              // Offset for each gene tree in the flattened array
    const int* orderingLengths,              // Length of each gene tree ordering
    int numCandidates,                       // Number of candidate bipartitions
    int numGeneTreeBips,                     // Number of gene tree bipartitions
    int maxTaxaCount                         // Maximum number of taxa
) {
    int candidateIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (candidateIdx >= numCandidates) return;
    
    double totalScore = 0.0;
    RangeBipartition candidate = candidates[candidateIdx];
    
    for (int i = 0; i < numGeneTreeBips; i++) {
        double score = calculateRangeScore(candidate, geneTreeBips[i],
                                          geneTreeOrderings, orderingOffsets,
                                          orderingLengths, maxTaxaCount);
        totalScore += score * frequencies[i];
    }
    
    weights[candidateIdx] = totalScore;
}

// Host function to launch the kernel
extern "C" {
    void launchRangeWeightCalculation(
        RangeBipartition* hCandidates,
        RangeBipartition* hGeneTreeBips,
        int* hFrequencies,
        double* hWeights,
        int** hGeneTreeOrderings,
        int* hOrderingLengths,
        int numCandidates,
        int numGeneTreeBips,
        int numGeneTrees,
        int maxTaxaCount
    ) {
        printf("Starting CUDA range weight calculation...\n");
        
        // Calculate total size needed for flattened orderings
        int totalOrderingSize = 0;
        for (int i = 0; i < numGeneTrees; i++) {
            totalOrderingSize += hOrderingLengths[i];
        }
        
        // Prepare flattened orderings and offsets
        int* flattenedOrderings = new int[totalOrderingSize];
        int* orderingOffsets = new int[numGeneTrees];
        
        int offset = 0;
        for (int i = 0; i < numGeneTrees; i++) {
            orderingOffsets[i] = offset;
            for (int j = 0; j < hOrderingLengths[i]; j++) {
                flattenedOrderings[offset + j] = hGeneTreeOrderings[i][j];
            }
            offset += hOrderingLengths[i];
        }
        
        // Device allocations
        RangeBipartition *dCandidates, *dGeneTreeBips;
        int *dFrequencies, *dOrderingLengths, *dOrderingOffsets;
        double *dWeights;
        int *dFlattenedOrderings;
        
        // Allocate device memory
        cudaMalloc(&dCandidates, numCandidates * sizeof(RangeBipartition));
        cudaMalloc(&dGeneTreeBips, numGeneTreeBips * sizeof(RangeBipartition));
        cudaMalloc(&dFrequencies, numGeneTreeBips * sizeof(int));
        cudaMalloc(&dWeights, numCandidates * sizeof(double));
        cudaMalloc(&dFlattenedOrderings, totalOrderingSize * sizeof(int));
        cudaMalloc(&dOrderingOffsets, numGeneTrees * sizeof(int));
        cudaMalloc(&dOrderingLengths, numGeneTrees * sizeof(int));
        
        // Copy data to device
        cudaMemcpy(dCandidates, hCandidates, numCandidates * sizeof(RangeBipartition), cudaMemcpyHostToDevice);
        cudaMemcpy(dGeneTreeBips, hGeneTreeBips, numGeneTreeBips * sizeof(RangeBipartition), cudaMemcpyHostToDevice);
        cudaMemcpy(dFrequencies, hFrequencies, numGeneTreeBips * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dFlattenedOrderings, flattenedOrderings, totalOrderingSize * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dOrderingOffsets, orderingOffsets, numGeneTrees * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dOrderingLengths, hOrderingLengths, numGeneTrees * sizeof(int), cudaMemcpyHostToDevice);
        
        // Launch kernel
        int blockSize = 256;
        int gridSize = (numCandidates + blockSize - 1) / blockSize;
        
        printf("Launching kernel with %d blocks of %d threads each\n", gridSize, blockSize);
        calculateRangeWeightsKernel<<<gridSize, blockSize>>>(
            dCandidates, dGeneTreeBips, dFrequencies, dWeights,
            dFlattenedOrderings, dOrderingOffsets, dOrderingLengths,
            numCandidates, numGeneTreeBips, maxTaxaCount
        );
        
        // Check for kernel launch errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA kernel launch error: %s\n", cudaGetErrorString(err));
        }
        
        // Wait for completion
        cudaDeviceSynchronize();
        
        // Copy results back to host
        cudaMemcpy(hWeights, dWeights, numCandidates * sizeof(double), cudaMemcpyDeviceToHost);
        
        // Cleanup device memory
        cudaFree(dCandidates);
        cudaFree(dGeneTreeBips);
        cudaFree(dFrequencies);
        cudaFree(dWeights);
        cudaFree(dFlattenedOrderings);
        cudaFree(dOrderingOffsets);
        cudaFree(dOrderingLengths);
        
        // Cleanup host memory
        delete[] flattenedOrderings;
        delete[] orderingOffsets;
        
        printf("CUDA range weight calculation completed\n");
    }
}
