#include <stdio.h>
#include <cuda_runtime.h>

#define N 16384     // Size of A
#define M 17408     // Size of B

// Example function f(a, b) = a & b
__device__ __forceinline__ int f(int a, int b) {
    return a & b;
}

// Kernel: each thread computes result[i] = sum_j f(A[j], B[i])
__global__ void compute_and_sums(int* A, int* B, int* result, int n, int m) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= m) return;

    int sum = 0;
    for (int j = 0; j < n; ++j) {
        sum += f(A[j], B[i]);
    }
    result[i] = sum;
}

int main() {
    // Host allocations
    int* hA = (int*) malloc(N * sizeof(int));
    int* hB = (int*) malloc(M * sizeof(int));
    int* hRes = (int*) malloc(M * sizeof(int));

    // Initialize with small integers to keep sums in safe range
    for (int i = 0; i < N; ++i) hA[i] = i % 64;
    for (int i = 0; i < M; ++i) hB[i] = (i * 3) % 128;

    // Device allocations
    int *dA, *dB, *dRes;
    cudaMalloc(&dA, N * sizeof(int));
    cudaMalloc(&dB, M * sizeof(int));
    cudaMalloc(&dRes, M * sizeof(int));

    // Copy data to GPU
    cudaMemcpy(dA, hA, N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dB, hB, M * sizeof(int), cudaMemcpyHostToDevice);

    // Launch kernel
    int blockSize = 256;
    int gridSize = (M + blockSize - 1) / blockSize;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    compute_and_sums<<<gridSize, blockSize>>>(dA, dB, dRes, N, M);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);

    // Copy result back
    cudaMemcpy(hRes, dRes, M * sizeof(int), cudaMemcpyDeviceToHost);

    // Print some results
    printf("First 5 results:\n");
    for (int i = 0; i < 5; ++i) {
        printf("Result[%d] = %d\n", i, hRes[i]);
    }
    printf("Time taken by GPU kernel: %.3f ms\n", ms);

    // Cleanup
    cudaFree(dA); cudaFree(dB); cudaFree(dRes);
    free(hA); free(hB); free(hRes);

    return 0;
}

