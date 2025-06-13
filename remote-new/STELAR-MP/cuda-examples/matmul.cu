#include <stdio.h>
#include <cuda_runtime.h>

#define N 512
#define TILE 16

__global__ void matmul_tiled(float* C, float* A, float* B) {
    __shared__ float tileA[TILE][TILE];
    __shared__ float tileB[TILE][TILE];

    int row = blockIdx.y * TILE + threadIdx.y;
    int col = blockIdx.x * TILE + threadIdx.x;
    float sum = 0.0f;

    for (int t = 0; t < N / TILE; ++t) {
        tileA[threadIdx.y][threadIdx.x] = A[row * N + t * TILE + threadIdx.x];
        tileB[threadIdx.y][threadIdx.x] = B[(t * TILE + threadIdx.y) * N + col];
        __syncthreads();

        for (int k = 0; k < TILE; ++k)
            sum += tileA[threadIdx.y][k] * tileB[k][threadIdx.x];

        __syncthreads();
    }

    C[row * N + col] = sum;
}

int main() {
    size_t size = N * N * sizeof(float);
    float *hA, *hB, *hC;
    cudaMallocHost(&hA, size);
    cudaMallocHost(&hB, size);
    cudaMallocHost(&hC, size);

    for (int i = 0; i < N * N; ++i) {
        hA[i] = 1.0f;
        hB[i] = 2.0f;
    }

    float *dA, *dB, *dC;
    cudaMalloc(&dA, size);
    cudaMalloc(&dB, size);
    cudaMalloc(&dC, size);

    cudaMemcpy(dA, hA, size, cudaMemcpyHostToDevice);
    cudaMemcpy(dB, hB, size, cudaMemcpyHostToDevice);

    dim3 block(TILE, TILE);
    dim3 grid(N / TILE, N / TILE);

    // === Benchmark Start ===
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    // Launch GPU kernel
    matmul_tiled<<<grid, block>>>(dC, dA, dB);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    // === Benchmark End ===

    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);

    cudaMemcpy(hC, dC, size, cudaMemcpyDeviceToHost);
    printf("Sample result hC[0] = %.2f\n", hC[0]);  // Should be 1024
    printf("Time taken by GPU kernel: %.3f ms\n", ms);

    cudaFree(dA); cudaFree(dB); cudaFree(dC);
    cudaFreeHost(hA); cudaFreeHost(hB); cudaFreeHost(hC);
    return 0;
}


