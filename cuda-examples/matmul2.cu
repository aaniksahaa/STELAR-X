#include <stdio.h>
#include <cuda_runtime.h>

#define N 512

__global__ void matmul_naive(float* C, float* A, float* B) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < N && col < N) {
        float sum = 0.0f;
        for (int k = 0; k < N; ++k) {
            sum += A[row * N + k] * B[k * N + col];
        }
        C[row * N + col] = sum;
    }
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

    dim3 block(16, 16);
    dim3 grid((N + 15) / 16, (N + 15) / 16);

    // === Benchmark Start ===
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    matmul_naive<<<grid, block>>>(dC, dA, dB);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    // === Benchmark End ===

    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);

    cudaMemcpy(hC, dC, size, cudaMemcpyDeviceToHost);
    printf("Sample result hC[0] = %.2f\n", hC[0]);  // Should be 1024
    printf("Time taken by naive GPU kernel: %.3f ms\n", ms);

    cudaFree(dA); cudaFree(dB); cudaFree(dC);
    cudaFreeHost(hA); cudaFreeHost(hB); cudaFreeHost(hC);
    return 0;
}

