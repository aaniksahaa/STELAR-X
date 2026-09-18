/**
 * stelarx_dp.cu — CUDA JNI for STELAR-X cross-tree DP transition discovery (Mode 2).
 *
 * Implements the size-binned hash-subtraction search:
 *   For every cluster A in X, for each size sz in [1..|A|/2],
 *   for each cluster B in bin(sz):
 *     residual = hash(A) - hash(B)    (arithmetic on raw sums/XORs)
 *     if residual ∈ X  →  emit transition  A → B | residual
 *
 * === GPU memory design (scalable) ===
 *
 * The cluster data (sums, XORs, sizes) for all |X| clusters is uploaded once
 * and stays on-device for the entire search.  Memory = O(|X| * m * 16) bytes.
 * For |X| = 1 M clusters with m=2 this is ~32 MB — fits in modern VRAM.
 * For even larger inputs (|X| >> 10 M), the cluster-data upload itself would
 * need to be batched; that extension is straightforward but not implemented here.
 *
 * The OUTPUT buffer is the hot-spot: the total number of found transitions can
 * be large (empirically O((nk)^1.73) in ASTRAL literature).  We avoid ever
 * allocating that entire buffer on the GPU at once by sub-batching:
 *
 *   For each size-sz round:
 *     sub-batch the "large clusters" (size >= 2*sz) into groups of at most
 *     (maxPerRound / binSzCount) clusters.  Each sub-batch launches one kernel
 *     whose output fits within a fixed GPU buffer of maxPerRound triples.
 *     After each sub-batch the buffer is flushed to the CPU and the GPU counter
 *     is reset to 0.
 *
 * GPU output buffer = maxPerRound * 3 * sizeof(int).
 *   Default maxPerRound = 10 000 000  →  120 MB.  Controllable from Java.
 *
 * === Hash table ===
 *
 * Open-addressing hash table with linear probing.
 * Load factor 0.25 (4× table size), max probe 512.
 * Slot layout: parallel arrays d_ht_sums[HT*m], d_ht_xors[HT*m], d_ht_sizes[HT],
 *              d_ht_idx[HT] (-1 = empty).
 * Slot hash = f(sums[0], size, m) using a 64-bit finalizer mix.
 * Memory = HT * (2*m*8 + 8) bytes.  For HT = 4*|X|, m=2: 160 bytes/cluster.
 * For |X| = 1 M: 160 MB.
 */

#include <jni.h>
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include "stelarx_platform.h"

// ─── Progress-bar helpers (host-side) ────────────────────────────────────────

static double dp_now_sec(void) {
    return stelarx_now_sec();
}

static void dp_fmt_duration(double secs, char* buf, int buflen) {
    int s = (int)secs;
    if (s < 60)
        snprintf(buf, buflen, "%ds", s);
    else if (s < 3600)
        snprintf(buf, buflen, "%dm%02ds", s / 60, s % 60);
    else
        snprintf(buf, buflen, "%dh%02dm", s / 3600, (s % 3600) / 60);
}

static int dp_use_color(void) {
    if (getenv("NO_COLOR"))    return 0;
    if (getenv("FORCE_COLOR")) return 1;
    return stelarx_stderr_isatty();
}

#define DP_BAR_W 28
static void dp_build_bar(char* buf, int done, int total) {
    int filled = (total > 0) ? (int)((double)done / total * DP_BAR_W + 0.5) : 0;
    if (filled > DP_BAR_W) filled = DP_BAR_W;
    int pos = 0;
    for (int i = 0; i < DP_BAR_W; i++) {
        if (i < filled) {
            buf[pos++] = '\xe2'; buf[pos++] = '\x96'; buf[pos++] = '\x88'; // █
        } else {
            buf[pos++] = '\xe2'; buf[pos++] = '\x96'; buf[pos++] = '\x91'; // ░
        }
    }
    buf[pos] = '\0';
}

// ─── Constants ───────────────────────────────────────────────────────────────
#define BLOCK_SIZE    256
#define MAX_M         4       // max supported hash seeds (m ≤ 4)
#define MAX_PROBE     512     // linear-probing limit; table at 0.25 load → rare
#define EMPTY_IDX     (-1)

// ─── Error checking ───────────────────────────────────────────────────────────
#define CUDA_CHECK(call, ret) do { \
    cudaError_t _e = (call); \
    if (_e != cudaSuccess) { \
        fprintf(stderr, "[stelarx_dp] CUDA error: %s  at %s:%d\n", \
                cudaGetErrorString(_e), __FILE__, __LINE__); \
        return (ret); \
    } \
} while(0)

// ─── Device: slot hash ────────────────────────────────────────────────────────
__device__ __forceinline__
unsigned int slotHash(const long long* sums, int size, int m, int HT_SIZE)
{
    // sums[0] is already pseudo-random (SplitMix64 sum).
    // Mix it with size and optionally sums[1] for extra entropy.
    unsigned long long h = (unsigned long long)sums[0];
    h ^= (unsigned long long)(unsigned int)size * 2654435761ULL;
    if (m > 1) h ^= (unsigned long long)sums[1] * 0x9e3779b97f4a7c15ULL;
    // Finalizer
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return (unsigned int)(h % (unsigned long long)(unsigned int)HT_SIZE);
}

// ─── Device: hash table match check ──────────────────────────────────────────
__device__ __forceinline__
bool htMatch(int s, int m, int qSize,
             const long long* qSums, const long long* qXors,
             const long long* ht_sums, const long long* ht_xors, const int* ht_sizes)
{
    if (ht_sizes[s] != qSize) return false;
    for (int seed = 0; seed < m; seed++) {
        if (ht_sums[s * m + seed] != qSums[seed]) return false;
        if (ht_xors[s * m + seed] != qXors[seed]) return false;
    }
    return true;
}

// ─── Device: hash table lookup → cluster index or EMPTY_IDX ──────────────────
__device__
int htLookup(const long long* qSums, const long long* qXors, int qSize, int m,
             const long long* ht_sums, const long long* ht_xors,
             const int* ht_sizes, const int* ht_idx, int HT_SIZE)
{
    unsigned int slot = slotHash(qSums, qSize, m, HT_SIZE);
    for (int probe = 0; probe < MAX_PROBE; probe++) {
        int s = (int)((slot + (unsigned int)probe) % (unsigned int)HT_SIZE);
        int idx = ht_idx[s];
        if (idx == EMPTY_IDX) return EMPTY_IDX;
        if (htMatch(s, m, qSize, qSums, qXors, ht_sums, ht_xors, ht_sizes))
            return idx;
    }
    return EMPTY_IDX;
}

// ─── Kernel: initialize hash table ───────────────────────────────────────────
__global__
void initHTKernel(int* ht_idx, int HT_SIZE)
{
    int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s < HT_SIZE) ht_idx[s] = EMPTY_IDX;
}

// ─── Kernel: insert all clusters into hash table ─────────────────────────────
// One thread per cluster c.
// Writes sums/xors/sizes BEFORE returning; ordering across threads is resolved
// by the cudaDeviceSynchronize() that separates build from query kernels.
__global__
void buildHTKernel(
    const long long* d_sums, const long long* d_xors, const int* d_sizes,
    int N, int m,
    long long* ht_sums, long long* ht_xors, int* ht_sizes, int* ht_idx, int HT_SIZE)
{
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= N) return;

    const long long* cSums = d_sums + (size_t)c * m;
    const long long* cXors = d_xors + (size_t)c * m;
    int cSize = d_sizes[c];

    unsigned int slot = slotHash(cSums, cSize, m, HT_SIZE);
    for (int probe = 0; probe < MAX_PROBE; probe++) {
        int s = (int)((slot + (unsigned int)probe) % (unsigned int)HT_SIZE);
        int old = atomicCAS(&ht_idx[s], EMPTY_IDX, c);
        if (old == EMPTY_IDX) {
            // Claimed slot s.  Write cluster data; the cudaDeviceSynchronize
            // after this kernel guarantees visibility to query kernels.
            for (int seed = 0; seed < m; seed++) {
                ht_sums[(size_t)s * m + seed] = cSums[seed];
                ht_xors[(size_t)s * m + seed] = cXors[seed];
            }
            ht_sizes[s] = cSize;
            return;
        }
    }
    // Table too full — should never happen at load factor 0.25
}

// ─── Kernel: cross-tree search for one sub-batch ─────────────────────────────
// Launch dimensions: aCount * binCount threads (1D).
// Each thread handles one (A, B) pair:
//   - computes residual = hash(A) − hash(B)
//   - looks up in hash table
//   - if found AND not a duplicate, atomically emits (idxA, idxB, idxRes)
//
// Duplicate avoidance (equal-size case):
//   When |B| == |residual| (only possible when |A| is even and |B| = |A|/2),
//   both (A, B, R) and (A, R, B) would be found.  Emit only when idxB < idxRes.
//
// Note: when |B| < |residual| (strict), A is never in the "large" set for
//   the round sz=|residual| (because |A| < 2*|residual| since |B| > 0),
//   so the symmetric pair is never processed — no extra check needed.
__global__
void crossTreeSearchKernel(
    const long long* d_sums, const long long* d_xors, const int* d_sizes,
    const int* d_aIdx, int aCount,
    const int* d_bIdx, int bCount,
    int sz, int m,
    const long long* ht_sums, const long long* ht_xors,
    const int* ht_sizes, const int* ht_idx, int HT_SIZE,
    int* d_outBuf, int* d_outCount, int maxOut)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= aCount * bCount) return;

    // Map 1D thread to (a, b) pair
    int ai = tid / bCount;
    int bi = tid % bCount;

    int idxA = d_aIdx[ai];
    int idxB = d_bIdx[bi];

    int szA = d_sizes[idxA];
    // szA >= 2*sz is guaranteed by how d_aIdx is selected (binStart[2*sz..N)),
    // but the check is cheap and prevents stray threads from doing wrong work.
    int resSize = szA - sz;
    if (resSize < sz) return;

    // Compute residual hash  (mod 2^64 arithmetic, C wraps naturally)
    long long rSums[MAX_M], rXors[MAX_M];
    for (int s = 0; s < m; s++) {
        rSums[s] = d_sums[(size_t)idxA * m + s] - d_sums[(size_t)idxB * m + s];
        rXors[s] = d_xors[(size_t)idxA * m + s] ^ d_xors[(size_t)idxB * m + s];
    }

    int idxRes = htLookup(rSums, rXors, resSize, m,
                          ht_sums, ht_xors, ht_sizes, ht_idx, HT_SIZE);
    if (idxRes == EMPTY_IDX) return;

    // Dedup equal-size pairs: emit only when idxB < idxRes
    if (resSize == sz && idxRes <= idxB) return;

    // Atomically reserve output slot
    int outPos = atomicAdd(d_outCount, 1);
    if (outPos < maxOut) {
        d_outBuf[outPos * 3 + 0] = idxA;
        d_outBuf[outPos * 3 + 1] = idxB;
        d_outBuf[outPos * 3 + 2] = idxRes;
    }
    // If outPos >= maxOut, the result is silently dropped.
    // The JNI caller sizes the buffer so this should not happen;
    // if it does, a warning is printed on the host side.
}

// ─── Kernel: reset outCount to 0 ─────────────────────────────────────────────
__global__ void resetCountKernel(int* d_outCount) {
    if (threadIdx.x == 0 && blockIdx.x == 0) *d_outCount = 0;
}

// ═══════════════════════════════════════════════════════════════════════════════
// JNI entry point
// ═══════════════════════════════════════════════════════════════════════════════
extern "C"
JNIEXPORT jintArray JNICALL
Java_stelarx_gpu_GPUDPBuilder_findCrossTreeTransitionsGPU(
    JNIEnv* env, jclass cls,
    jlongArray jClusterSums,        // [N*m]
    jlongArray jClusterXors,        // [N*m]
    jintArray  jClusterSizes,       // [N]
    jint       jN,
    jint       jM,
    jintArray  jSortedBySize,       // [N]
    jintArray  jBinStart,           // [maxSize+2]
    jint       jMaxSize,
    jint       jMaxPerRound,        // GPU output buffer size per sub-batch
    jdouble    jProgressInterval,   // min seconds between progress bar updates
    jint       jProgressMaxSteps)   // max number of progress bar lines
{
    const int    N              = (int)jN;
    const int    mSeeds         = (int)jM;
    const int    maxSz          = (int)jMaxSize;
    const int    maxPerRound    = (int)jMaxPerRound;
    const double progressInterval  = (double)jProgressInterval;
    const int    progressMaxSteps  = (int)jProgressMaxSteps;
    const int HT_SIZE    = (int)(std::max(16, N * 4)); // load factor 0.25

    if (N == 0) {
        jintArray empty = env->NewIntArray(1);
        int zero = 0;
        env->SetIntArrayRegion(empty, 0, 1, &zero);
        return empty;
    }

    // ── Pin host arrays ───────────────────────────────────────────────────────
    jlong* hSums     = env->GetLongArrayElements(jClusterSums,  NULL);
    jlong* hXors     = env->GetLongArrayElements(jClusterXors,  NULL);
    jint*  hSizes    = env->GetIntArrayElements( jClusterSizes, NULL);
    jint*  hSorted   = env->GetIntArrayElements( jSortedBySize, NULL);
    jint*  hBinStart = env->GetIntArrayElements( jBinStart,     NULL);

#define RELEASE_AND_RETURN(ret) do { \
    env->ReleaseLongArrayElements(jClusterSums,  hSums,     JNI_ABORT); \
    env->ReleaseLongArrayElements(jClusterXors,  hXors,     JNI_ABORT); \
    env->ReleaseIntArrayElements( jClusterSizes, hSizes,    JNI_ABORT); \
    env->ReleaseIntArrayElements( jSortedBySize, hSorted,   JNI_ABORT); \
    env->ReleaseIntArrayElements( jBinStart,     hBinStart, JNI_ABORT); \
    return (ret); \
} while(0)

    // ── Device allocations: cluster data ─────────────────────────────────────
    long long *d_sums = NULL, *d_xors = NULL;
    int *d_sizes = NULL;
    if (cudaMalloc(&d_sums,  (size_t)N * mSeeds * sizeof(long long)) != cudaSuccess ||
        cudaMalloc(&d_xors,  (size_t)N * mSeeds * sizeof(long long)) != cudaSuccess ||
        cudaMalloc(&d_sizes, (size_t)N * sizeof(int))                != cudaSuccess)
    {
        fprintf(stderr, "[stelarx_dp] cudaMalloc failed for cluster data\n");
        cudaFree(d_sums); cudaFree(d_xors); cudaFree(d_sizes);
        RELEASE_AND_RETURN(NULL);
    }
    cudaMemcpy(d_sums,  hSums,  (size_t)N * mSeeds * sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(d_xors,  hXors,  (size_t)N * mSeeds * sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sizes, hSizes, (size_t)N * sizeof(int),                cudaMemcpyHostToDevice);

    // ── Device allocations: hash table ───────────────────────────────────────
    long long *d_ht_sums = NULL, *d_ht_xors = NULL;
    int *d_ht_sizes = NULL, *d_ht_idx = NULL;
    if (cudaMalloc(&d_ht_sums,  (size_t)HT_SIZE * mSeeds * sizeof(long long)) != cudaSuccess ||
        cudaMalloc(&d_ht_xors,  (size_t)HT_SIZE * mSeeds * sizeof(long long)) != cudaSuccess ||
        cudaMalloc(&d_ht_sizes, (size_t)HT_SIZE * sizeof(int))                != cudaSuccess ||
        cudaMalloc(&d_ht_idx,   (size_t)HT_SIZE * sizeof(int))                != cudaSuccess)
    {
        fprintf(stderr, "[stelarx_dp] cudaMalloc failed for hash table\n");
        cudaFree(d_sums); cudaFree(d_xors); cudaFree(d_sizes);
        cudaFree(d_ht_sums); cudaFree(d_ht_xors); cudaFree(d_ht_sizes); cudaFree(d_ht_idx);
        RELEASE_AND_RETURN(NULL);
    }

    // Init + build hash table
    initHTKernel<<<(HT_SIZE + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(d_ht_idx, HT_SIZE);
    cudaDeviceSynchronize();
    buildHTKernel<<<(N + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(
        d_sums, d_xors, d_sizes, N, mSeeds,
        d_ht_sums, d_ht_xors, d_ht_sizes, d_ht_idx, HT_SIZE);
    cudaDeviceSynchronize();

    // ── Device allocations: per-round index arrays and output buffer ──────────
    int *d_aIdx = NULL, *d_bIdx = NULL;
    int *d_outBuf = NULL, *d_outCount = NULL;
    if (cudaMalloc(&d_aIdx,     (size_t)N * sizeof(int))          != cudaSuccess ||
        cudaMalloc(&d_bIdx,     (size_t)N * sizeof(int))          != cudaSuccess ||
        cudaMalloc(&d_outBuf,   (size_t)maxPerRound * 3 * sizeof(int)) != cudaSuccess ||
        cudaMalloc(&d_outCount, sizeof(int))                      != cudaSuccess)
    {
        fprintf(stderr, "[stelarx_dp] cudaMalloc failed for per-round buffers\n");
        cudaFree(d_sums); cudaFree(d_xors); cudaFree(d_sizes);
        cudaFree(d_ht_sums); cudaFree(d_ht_xors); cudaFree(d_ht_sizes); cudaFree(d_ht_idx);
        cudaFree(d_aIdx); cudaFree(d_bIdx); cudaFree(d_outBuf); cudaFree(d_outCount);
        RELEASE_AND_RETURN(NULL);
    }
    {   // zero the output counter
        int zero = 0;
        cudaMemcpy(d_outCount, &zero, sizeof(int), cudaMemcpyHostToDevice);
    }

    // Analytical VRAM budget: log all allocations for this phase
    {
        size_t clusterDataBytes = (size_t)N * mSeeds * 2 * sizeof(long long)
                                + (size_t)N * sizeof(int);
        size_t htBytes = (size_t)HT_SIZE * mSeeds * 2 * sizeof(long long)
                       + (size_t)HT_SIZE * 2 * sizeof(int);
        size_t outBufBytes = (size_t)maxPerRound * 3 * sizeof(int);
        size_t idxBufBytes = (size_t)N * 2 * sizeof(int);
        size_t totalBytes  = clusterDataBytes + htBytes + outBufBytes + idxBufBytes;
        size_t freeNow = 0, totalVRAM = 0;
        cudaMemGetInfo(&freeNow, &totalVRAM);
        fprintf(stderr,
            "[STELAR-X GPU] DP cross-tree allocations:\n"
            "  N=%d  HT_SIZE=%d  mSeeds=%d  maxPerRound=%d\n"
            "  cluster data : %6.1f MB\n"
            "  hash table   : %6.1f MB\n"
            "  output buf   : %6.1f MB  (%d triples max)\n"
            "  index bufs   : %6.1f MB\n"
            "  ─────────────────────────────────────────\n"
            "  DP total     : %6.1f MB   (VRAM free now: %.1f MB / %.1f MB)\n",
            N, HT_SIZE, mSeeds, maxPerRound,
            clusterDataBytes / 1e6,
            htBytes          / 1e6,
            outBufBytes      / 1e6, maxPerRound,
            idxBufBytes      / 1e6,
            totalBytes       / 1e6,
            freeNow / 1e6, totalVRAM / 1e6);
        fflush(stderr);
    }

    // ── Host accumulator for all found triples ────────────────────────────────
    std::vector<int> accum; // (idxA, idxB, idxRes) interleaved
    accum.reserve(std::min(N * 10, 1000000));

    // ── Count active size bins for progress reporting ─────────────────────────
    int activeBins = 0;
    for (int sz = 1; sz * 2 <= maxSz; sz++) {
        int bCount      = hBinStart[sz + 1] - hBinStart[sz];
        int twoSz       = 2 * sz;
        int aStart      = (twoSz <= maxSz + 1) ? hBinStart[twoSz] : N;
        int aTotalCount = N - aStart;
        if (bCount > 0 && aTotalCount > 0) activeBins++;
    }

    fprintf(stderr,
        "[STELAR-X GPU] DP cross-tree search starting:\n"
        "  clusters     : %d\n"
        "  active bins  : %d  (size range 1..%d)\n"
        "  output buf   : %d triples / sub-batch  (cap = %.0f MB)\n",
        N, activeBins, maxSz / 2,
        maxPerRound, (double)maxPerRound * 3 * sizeof(int) / 1e6);
    fflush(stderr);

    // ── Round-by-round size-binned search ─────────────────────────────────────
    // Two mutually exclusive progress modes (chosen at runtime):
    //   progressMaxSteps > 0  →  step mode:  print every (100/maxSteps)% advancement
    //   progressMaxSteps == 0 →  time mode:  print every progressInterval seconds
    const bool step_mode = (progressMaxSteps > 0);
    const char* GRN = dp_use_color() ? "\033[32m" : "";
    const char* RST = dp_use_color() ? "\033[0m"  : "";
    int    activeDone       = 0;
    double t_loop_start     = dp_now_sec();
    double t_last_print     = t_loop_start - progressInterval; // force first eligible (time mode)
    double last_pct_printed = step_mode ? -(100.0 / progressMaxSteps) : 0.0; // force first eligible (step mode)
    char   dp_bar_buf[DP_BAR_W * 3 + 1];

    for (int sz = 1; sz * 2 <= maxSz; sz++) {
        // bin[sz]: sortedBySize[ binStart[sz] .. binStart[sz+1] )
        int bStart = hBinStart[sz];
        int bEnd   = hBinStart[sz + 1];
        int bCount = bEnd - bStart;
        if (bCount == 0) continue;

        // large clusters: size >= 2*sz → sortedBySize[ binStart[2*sz] .. N )
        int twoSz   = 2 * sz;
        int aStart  = (twoSz <= maxSz + 1) ? (int)hBinStart[twoSz] : N;
        int aTotalCount = N - aStart;
        if (aTotalCount == 0) continue;

        activeDone++;

        // Upload bin[sz] indices to d_bIdx
        cudaMemcpy(d_bIdx, hSorted + bStart, (size_t)bCount * sizeof(int),
                   cudaMemcpyHostToDevice);

        // Sub-batch the A clusters so each kernel's output fits in maxPerRound
        int batchA = (bCount > 0) ? std::max(1, maxPerRound / bCount) : aTotalCount;

        for (int aOff = 0; aOff < aTotalCount; aOff += batchA) {
            int aCount = std::min(batchA, aTotalCount - aOff);

            // Upload this sub-batch of A indices
            cudaMemcpy(d_aIdx, hSorted + aStart + aOff,
                       (size_t)aCount * sizeof(int), cudaMemcpyHostToDevice);

            // Reset output counter
            resetCountKernel<<<1, 1>>>(d_outCount);
            cudaDeviceSynchronize();

            // Launch search kernel
            int totalThreads = aCount * bCount;
            int numBlocks    = (totalThreads + BLOCK_SIZE - 1) / BLOCK_SIZE;
            crossTreeSearchKernel<<<numBlocks, BLOCK_SIZE>>>(
                d_sums, d_xors, d_sizes,
                d_aIdx, aCount,
                d_bIdx, bCount,
                sz, mSeeds,
                d_ht_sums, d_ht_xors, d_ht_sizes, d_ht_idx, HT_SIZE,
                d_outBuf, d_outCount, maxPerRound);
            cudaDeviceSynchronize();

            // Copy back result count
            int outCount = 0;
            cudaMemcpy(&outCount, d_outCount, sizeof(int), cudaMemcpyDeviceToHost);

            if (outCount > maxPerRound) {
                int lost = outCount - maxPerRound;
                fprintf(stderr,
                    "\n"
                    "  ╔══════════════════════════════════════════════════════════════╗\n"
                    "  ║           !!!  CRITICAL WARNING  !!!                        ║\n"
                    "  ╠══════════════════════════════════════════════════════════════╣\n"
                    "  ║  GPU OUTPUT BUFFER OVERFLOW IN CROSS-TREE DP PHASE          ║\n"
                    "  ║                                                              ║\n"
                    "  ║  Sub-batch produced %d triples but the GPU output\n"
                    "  ║  buffer cap is %d (%.0f MB).  %d transitions\n"
                    "  ║  have been SILENTLY DROPPED.                                 ║\n"
                    "  ║                                                              ║\n"
                    "  ║  CONSEQUENCE: the DP search space is INCOMPLETE.  Missing   ║\n"
                    "  ║  candidate splits may degrade the inferred species tree.     ║\n"
                    "  ║                                                              ║\n"
                    "  ║  Context:  size-bin sz=%d,  aOff=%d                          ║\n"
                    "  ║                                                              ║\n"
                    "  ║  TO FIX: rerun with --gpu-dp-state-space-construction-output-cap  ║\n"
                    "  ║  set to a larger value, e.g.  --gpu-dp-state-space-...  1g  ║\n"
                    "  ║  or report this dataset to the STELAR-X maintainers.        ║\n"
                    "  ╚══════════════════════════════════════════════════════════════╝\n"
                    "\n",
                    outCount,
                    maxPerRound, (double)maxPerRound * 3 * sizeof(int) / 1e6, lost,
                    sz, aOff);
                fflush(stderr);
                outCount = maxPerRound;
            }

            if (outCount > 0) {
                int prev = (int)accum.size();
                accum.resize(prev + outCount * 3);
                cudaMemcpy(accum.data() + prev, d_outBuf,
                           (size_t)outCount * 3 * sizeof(int), cudaMemcpyDeviceToHost);
            }
        }

        // ── Progress bar ──────────────────────────────────────────────────────
        if (activeBins > 1) {
            int    rem     = activeBins - activeDone;
            bool   is_last = (rem == 0);
            double pct     = 100.0 * activeDone / activeBins;
            double now     = dp_now_sec();
            bool   should_print;
            if (step_mode)
                should_print = is_last || (pct - last_pct_printed >= 100.0 / progressMaxSteps);
            else
                should_print = is_last || (now - t_last_print >= progressInterval);
            if (should_print) {
                t_last_print     = now;
                last_pct_printed = pct;
                double elapsed   = now - t_loop_start;
                int    found     = (int)(accum.size() / 3);
                dp_build_bar(dp_bar_buf, activeDone, activeBins);
                char dur_buf[32];
                dp_fmt_duration(elapsed, dur_buf, sizeof(dur_buf));
                if (is_last) {
                    // Final: end with newline, pad to clear any previous line
                    fprintf(stderr,
                        "\r  %s[GPU]%s dp      %s[%s]%s  %d/%d  100%%  %-8s  found: %-12d\n",
                        GRN, RST, GRN, dp_bar_buf, RST,
                        activeBins, activeBins, dur_buf, found);
                } else {
                    // Intermediate: \r with fixed-width fields so line never shrinks
                    fprintf(stderr,
                        "\r  %s[GPU]%s dp      %s[%s]%s  %4d/%-4d  %5.1f%%  elapsed: %-8s  found: %-12d",
                        GRN, RST, GRN, dp_bar_buf, RST,
                        activeDone, activeBins, pct, dur_buf, found);
                }
                fflush(stderr);
            }
        }
    }

    // Final one-liner when only a single active bin (no bar was printed)
    if (activeBins <= 1) {
        double elapsed = dp_now_sec() - t_loop_start;
        char dur_buf[32];
        dp_fmt_duration(elapsed, dur_buf, sizeof(dur_buf));
        fprintf(stderr,
            "  %s[GPU]%s dp      done in %s  found: %d transitions\n",
            GRN, RST, dur_buf, (int)(accum.size() / 3));
        fflush(stderr);
    }

    // ── Cleanup device memory ─────────────────────────────────────────────────
    cudaFree(d_sums);     cudaFree(d_xors);     cudaFree(d_sizes);
    cudaFree(d_ht_sums);  cudaFree(d_ht_xors);  cudaFree(d_ht_sizes); cudaFree(d_ht_idx);
    cudaFree(d_aIdx);     cudaFree(d_bIdx);
    cudaFree(d_outBuf);   cudaFree(d_outCount);

    // ── Build jintArray result: [totalCount, triples...] ─────────────────────
    int totalPairs = (int)(accum.size() / 3);
    jintArray result = env->NewIntArray(1 + totalPairs * 3);
    if (result != NULL) {
        env->SetIntArrayRegion(result, 0, 1, &totalPairs);
        if (totalPairs > 0)
            env->SetIntArrayRegion(result, 1, totalPairs * 3, accum.data());
    }

    RELEASE_AND_RETURN(result);
#undef RELEASE_AND_RETURN
}
